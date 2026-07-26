#!/usr/bin/env python3
"""Deterministic DEM fixtures with a depression pit placed relative to a tile seam.

These exercise the one case Barnes' halo-free perimeter-strip join may miss for the
depression *hierarchy* (PARALLEL_DEPHIER_PLAN.md §6.6): a basin whose lowest cells
straddle a tile boundary. When the domain is split at `--split-col`, both halves of
such a bowl drain toward the seam, both become BOUNDARY, and the basin can appear in
neither tile's local depression list. If the stitch reproduces the serial depression
there, the strip suffices; if the oracle shows a missing depression, we need a halo.

Modes (bowl x-centre relative to a split between columns split_col-1 and split_col):
  seam      -- centre on the seam line (x = split_col - 0.5): the two lowest cells are
               split_col-1 and split_col, one in each tile. THE STRESS CASE.
  offset    -- centre one cell inside the left tile (x = split_col - 1): pit adjacent
               to, but not on, the seam.
  interior  -- centre well inside the left tile (x = split_col // 2): CONTROL; the
               left tile forms the depression normally.

The bowl is a cone dug into a flat plateau inside a one-cell ocean ring (NODATA), so
there is exactly one interior depression whose pit is at the chosen location and whose
outlet is the plateau rim. Emits the ESRI ASCII grid dephier reads natively, and
prints a machine-readable FIXTURE line with the pit location and split column.
"""
import argparse
import sys

import numpy as np


def bowl_dem(size, split_col, mode, plateau=100.0, pit_elev=1.0, slope=3.0, nodata=-9999.0):
    if mode == "seam":
        cx = split_col - 0.5
    elif mode == "offset":
        cx = float(split_col - 1)
    elif mode == "interior":
        cx = float(split_col // 2)
    else:
        sys.exit(f"unknown mode {mode!r}")
    cy = size // 2

    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    dist = np.hypot(xx - cx, yy - cy)
    # Cone rising from the pit at `slope` m/cell, clipped to the surrounding plateau.
    topo = np.minimum(plateau, pit_elev + slope * dist)

    # One-cell ocean ring (drainage boundary) as NODATA.
    topo[0, :] = topo[-1, :] = topo[:, 0] = topo[:, -1] = nodata
    return topo, cx, cy


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
    ap.add_argument("--size", type=int, default=41, help="grid edge length (cells); square")
    ap.add_argument("--split-col", type=int, default=None,
                    help="column of the future tile seam (default: size//2)")
    ap.add_argument("--mode", choices=["seam", "offset", "interior"], default="seam")
    ap.add_argument("-o", "--out", required=True, help="output .dem path")
    a = ap.parse_args()

    split_col = a.split_col if a.split_col is not None else a.size // 2
    nodata = -9999.0
    topo, cx, cy = bowl_dem(a.size, split_col, a.mode, nodata=nodata)

    # Verify the fixture actually places the minimum land cell where intended.
    land = np.where(topo == nodata, np.inf, topo)
    min_val = land.min()
    min_cells = np.argwhere(land == min_val)  # (row, col) pairs
    min_cols = sorted({int(c) for _, c in min_cells})
    expected = {
        "seam":     [split_col - 1, split_col],
        "offset":   [split_col - 1],
        "interior": [split_col // 2],
    }[a.mode]
    if min_cols != expected:
        sys.exit(f"FIXTURE CHECK FAILED: min at cols {min_cols}, expected {expected} "
                 f"(mode={a.mode}, split_col={split_col})")

    write_esri_ascii(a.out, topo, nodata=nodata)
    print(f"FIXTURE {a.out} size={a.size} split_col={split_col} mode={a.mode} "
          f"pit_cols={min_cols} pit_row={cy} pit_elev={min_val:g}")


if __name__ == "__main__":
    main()
