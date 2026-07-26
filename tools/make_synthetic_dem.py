#!/usr/bin/env python3
"""Generate a synthetic DEM (ESRI ASCII grid) for DepressionHierarchy footprint sweeps.

Reuses WTM/tests/spectral_terrain.fractal_terrain -- the 1/|k|^beta spectral
generator -- as the section-8 precursor terrain source (PARALLEL_DEPHIER_PLAN.md).

IMPORTANT (the "pleasing number" caveat from section 8): for footprint *bounding*
we want the ADVERSARIAL, count-maximizing regime -- low beta, no kmax (rough) --
NOT the smooth, band-limited posture the generator defaults to for FSM golden
tests. Band-limiting suppresses small depressions and would understate the count.

Ocean encoding: the generator's one-cell ocean ring (mask==0) is written as
NODATA, which GetDepressionHierarchy treats as OCEAN (the drainage boundary).

The WTM checkout location defaults to /home/awickert/models/WTM; override with
the WTM_DIR environment variable.
"""
import argparse
import os
import sys

import numpy as np


def load_fractal_terrain():
    wtm = os.environ.get("WTM_DIR", "/home/awickert/models/WTM")
    tests = os.path.join(wtm, "tests")
    if not os.path.isfile(os.path.join(tests, "spectral_terrain.py")):
        sys.exit(f"spectral_terrain.py not found under {tests!r}; set WTM_DIR to your WTM checkout.")
    sys.path.insert(0, tests)
    from spectral_terrain import fractal_terrain
    return fractal_terrain


def write_esri_ascii(path, topo, nodata, cellsize=1.0):
    """Write a 2D float array as an ESRI ASCII grid (the test_cases/*.dem format)."""
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
    ap.add_argument("--size", type=int, default=512, help="grid edge length (cells); square")
    ap.add_argument("--beta", type=float, default=1.5,
                    help="spectral exponent; LOWER = rougher = more depressions (default 1.5)")
    ap.add_argument("--relief", type=float, default=1000.0, help="peak-to-peak relief (m)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed")
    ap.add_argument("--kmax", type=float, default=None,
                    help="band-limit wavenumber; leave UNSET for the adversarial count (see module docstring)")
    ap.add_argument("-o", "--out", required=True, help="output .dem path")
    a = ap.parse_args()

    fractal_terrain = load_fractal_terrain()

    topo, mask = fractal_terrain(
        (a.size, a.size), beta=a.beta, relief=a.relief,
        base=1.0, seed=a.seed, kmax=a.kmax, ocean_ring=True,
    )

    nodata = -9999.0
    topo = topo.astype(np.float64)
    topo[mask == 0] = nodata

    write_esri_ascii(a.out, topo, nodata=nodata)

    land = topo[mask != 0]
    print(f"wrote {a.out}: {a.size}x{a.size}, beta={a.beta}, seed={a.seed}, "
          f"kmax={a.kmax}, land={land.size}, elev[{land.min():.1f},{land.max():.1f}]")


if __name__ == "__main__":
    main()
