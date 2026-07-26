#!/usr/bin/env python3
"""Deterministic DEM fixtures with depressions whose pits are safely INTERIOR to
their tiles -- the benign case to validate the stitch first, before tackling
pits-on-edges (PARALLEL_DEPHIER_PLAN.md §6.6, and see make_edge_depression_dem.py).

Terrain = a gently TILTED base plane (so every cell elevation is distinct: no ties,
hence a unique depression tree and a clean serial-vs-stitch bit-identity check) inside
a one-cell ocean ring (NODATA), with cone bowls dug at interior locations chosen to
stay clear of both the ocean ring and the future tile seam.

Scenarios (domain split between columns split_col-1 and split_col):
  single    -- one bowl well inside the left tile.
  two       -- one bowl inside each tile; exercises namespace remap + reassembly of
               WHOLE (unsplit) depressions, both spilling independently to the ocean.
  connected -- one bowl inside each tile, joined by a trench whose sill sits ON the
               seam and below the ocean-spill, so the two merge into a meta-depression
               across the seam. Pits stay interior; only the OUTLET crosses the seam --
               the key cross-tile-merge stitch test, with no pit on an edge.

Self-check: computes the actual local minima of the generated DEM and asserts each
is a bowl centre, is off the ocean ring, and is at least `edge_margin` cells from the
seam -- i.e. genuinely no pit on an edge. Prints a FIXTURE line with pit locations.
"""
import argparse
import sys

import numpy as np


def build(size, split_col, scenario, edge_margin=3,
          base=25.0, gx=0.02, gy=0.03, slope=4.0, pit_elev=1.0, nodata=-9999.0):
    row = size // 2
    if scenario == "single":
        centres = [(row, split_col // 2)]
    elif scenario in ("two", "connected"):
        centres = [(row, split_col // 2),
                   (row, split_col + (size - split_col) // 2)]
    else:
        sys.exit(f"unknown scenario {scenario!r}")

    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    topo = base + gx * xx + gy * yy               # tilted plane: distinct, no ties
    for (cy, cx) in centres:
        dist = np.hypot(xx - cx, yy - cy)
        topo = np.minimum(topo, pit_elev + slope * dist)  # carve a cone bowl

    if scenario == "connected":
        # A trench along the pit row joining the two bowls, with its high point (the
        # inter-bowl saddle) on the seam and below the ocean-spill. `sill` sets the
        # saddle height; the small +tilt keeps every cell distinct (no ties).
        cxl, cxr = centres[0][1], centres[1][1]
        sill, trench_floor = 12.0, 6.0
        for x in range(cxl, cxr + 1):
            ramp = trench_floor + (sill - trench_floor) * (1.0 - abs(x - split_col) / (cxr - cxl))
            for y in (row - 1, row, row + 1):
                topo[y, x] = min(topo[y, x], ramp + gx * x + gy * y)

    topo[0, :] = topo[-1, :] = topo[:, 0] = topo[:, -1] = nodata
    return topo, centres


def local_minima(topo, nodata):
    """(row, col) of depression pits: land cells strictly lower than all 8 land
    neighbours AND not adjacent to the ocean ring. A cell touching NODATA drains to
    the ocean (it is coastal), so it is not a depression pit even if it is a local
    low of the tilted plane."""
    ny, nx = topo.shape
    mins = []
    for y in range(1, ny - 1):
        for x in range(1, nx - 1):
            v = topo[y, x]
            if v == nodata:
                continue
            lowest = True
            touches_ocean = False
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dy == 0 and dx == 0:
                        continue
                    n = topo[y + dy, x + dx]
                    if n == nodata:
                        touches_ocean = True
                    elif n < v:
                        lowest = False
            if lowest and not touches_ocean:
                mins.append((y, x))
    return mins


def write_esri_ascii(path, topo, nodata, cellsize=1.0):
    ny, nx = topo.shape
    with open(path, "w") as f:
        f.write(f"ncols         {nx}\n")
        f.write(f"nrows         {ny}\n")
        f.write("xllcorner     0\n")
        f.write("yllcorner     0\n")
        f.write(f"cellsize      {cellsize}\n")
        f.write(f"NODATA_value  {nodata:g}\n")
        np.savetxt(f, topo, fmt="%.6g")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--size", type=int, default=41)
    ap.add_argument("--split-col", type=int, default=None, help="future tile seam (default size//2)")
    ap.add_argument("--scenario", choices=["single", "two", "connected"], default="two")
    ap.add_argument("--edge-margin", type=int, default=3, help="min cells a pit must be from the seam")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    split_col = a.split_col if a.split_col is not None else a.size // 2
    nodata = -9999.0
    topo, centres = build(a.size, split_col, a.scenario, edge_margin=a.edge_margin, nodata=nodata)

    mins = local_minima(topo, nodata)
    # Every local minimum must be a requested bowl centre (no spurious pits from the plane).
    if sorted(mins) != sorted(centres):
        sys.exit(f"FIXTURE CHECK FAILED: local minima {sorted(mins)} != bowl centres {sorted(centres)}")
    # No pit on/near the seam or ring.
    for (r, c) in mins:
        if abs(c - split_col) < a.edge_margin or abs(c - (split_col - 1)) < a.edge_margin:
            sys.exit(f"FIXTURE CHECK FAILED: pit at col {c} within {a.edge_margin} of seam {split_col}")
        if r in (1, a.size - 2) or c in (1, a.size - 2):
            sys.exit(f"FIXTURE CHECK FAILED: pit at ({r},{c}) touches the ocean ring")

    write_esri_ascii(a.out, topo, nodata=nodata)
    pits = ";".join(f"{r},{c}" for (r, c) in sorted(mins))
    print(f"FIXTURE {a.out} size={a.size} split_col={split_col} scenario={a.scenario} "
          f"npits={len(mins)} pits(row,col)={pits}")


if __name__ == "__main__":
    main()
