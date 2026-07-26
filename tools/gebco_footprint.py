#!/usr/bin/env python3
"""Measure DepressionHierarchy footprint on real GEBCO 30" tiles (PARALLEL_DEPHIER_PLAN.md section 8).

The distributed build distributes the O(N) flood (Phase A+B) but *centralizes* the
light tree (Phase C) on one node. That is only viable if #depressions + boundary << N
on real terrain. This tool measures depression *density* across contrasting terrains so
the centralized-Phase-C decision rests on data, not on smooth synthetic fractals.
See PARALLEL_DEPHIER_ENGINEERING.md section 4 for the recorded verdict.

Pipeline per tile: gdal_translate -srcwin (cut from the global GeoTIFF) -> preprocess
(sea z<=0 and the 1-cell border ring -> ocean, so an inland tile has a base level) ->
dephier_stats (unmodified serial GetDepressionHierarchy; emits a DHSTATS line).

Prereqs:
  - GEBCO 30" as GeoTIFF. The distributed archive is legacy GMT 1-D netCDF; convert once:
        gmt grdconvert gebco_08.nc gebco_08.tif=gd:GTiff
    Grid: 43200x21600, origin (-180,90), 30" (1/120 deg) pixels.
  - dephier_stats.exe built (cmake -DUSE_GDAL=ON).
  - gdal_translate on PATH; numpy. (osgeo not required -- we route through AAIGrid.)
  - Run under the libstdc++/GDAL LD path, e.g.:
        LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$CONDA_PREFIX/lib

Usage:
    python3 gebco_footprint.py <gebco.tif> <dephier_stats.exe> [out.txt]

srcwin from a lon/lat box (whole-globe GeoTIFF, row 0 = 90 N):
    xoff = (lon_min + 180) * 120     xsize = (lon_max - lon_min) * 120
    yoff = (90 - lat_max) * 120      ysize = (lat_max - lat_min) * 120
"""
import subprocess, sys, os, tempfile
import numpy as np

# name -> srcwin (xoff, yoff, xsize, ysize), contrasting terrains; 15deg-wide unless noted
REGIONS = {
    "andes":            (12600, 13200, 1800, 1800),  # high relief
    "tibet":            (31440,  7440, 1800, 1200),  # high relief, endorheic plateau
    "shield":           ( 9600,  3600, 1800, 1200),  # Canadian Shield, lake-dense
    "sahara":           (22560,  7440, 1800, 1200),  # desert flats
    "amazon":           (13200, 10560, 1800, 1200),  # low relief, wet
    "centralus":        ( 9000,  5400, 1800, 1200),
    "fennoscandia":     (23040,  2640, 1800, 1200),  # glacial lakes
    "australia":        (36960, 12960, 1800, 1200),  # interior, endorheic
    "wsiberia":         (28800,  3000, 1800, 1200),  # flat/wet
    "greatbasin":       ( 7320,  5760, 1080,  720),  # endorheic, ~9deg
    "centralasia_big":  (27600,  4800, 3600, 2400),  # 30x20deg, tile-size scaling check
}
ND = -9999.0

def preprocess(asc_in, asc_out):
    Z = np.loadtxt(asc_in, skiprows=6).astype(np.float64)
    if Z.ndim == 1:
        Z = Z.reshape(1, -1)
    Z[np.isnan(Z)] = ND
    Z[Z <= 0.0] = ND                                  # sea / at-or-below sea level -> ocean
    Z[0, :] = Z[-1, :] = Z[:, 0] = Z[:, -1] = ND      # open border ring -> base level
    nr, nc = Z.shape
    with open(asc_out, "w") as f:
        f.write(f"ncols {nc}\nnrows {nr}\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value {int(ND)}\n")
        np.savetxt(f, Z, fmt="%.1f")
    return int(np.sum(Z != ND))

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    gebco, stats_exe = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else "gebco_footprint.txt"
    with open(out, "w") as log, tempfile.TemporaryDirectory() as td:
        tot_dep = tot_land = 0
        for name, (xo, yo, xs, ys) in REGIONS.items():
            raw = os.path.join(td, f"{name}.asc"); pp = os.path.join(td, f"{name}_pp.asc")
            if subprocess.run(["gdal_translate", "-q", "-of", "AAIGrid", "-srcwin",
                               str(xo), str(yo), str(xs), str(ys), gebco, raw]).returncode:
                log.write(f"{name} EXTRACT_FAIL\n"); continue
            land = preprocess(raw, pp)
            r = subprocess.run([stats_exe, pp, "-99999"], capture_output=True, text=True)
            line = next((l for l in r.stdout.splitlines() if l.startswith("DHSTATS")), None)
            if not line:
                log.write(f"{name} SKIPPED (no ocean / error)\n"); continue
            ndep = int(line.split()[6])
            tot_dep += ndep; tot_land += land
            log.write(f"{name} density={ndep/land:.4f} {line}\n"); log.flush()
            print(f"{name:16s} land={land:>9d} dep={ndep:>8d} density={ndep/land:6.2%}")
        if tot_land:
            print(f"\naggregate density = {tot_dep/tot_land:.2%} ({tot_dep} dep / {tot_land} land cells)")
            log.write(f"AGGREGATE density={tot_dep/tot_land:.4f} dep={tot_dep} land={tot_land}\n")

if __name__ == "__main__":
    main()
