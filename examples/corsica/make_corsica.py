#!/usr/bin/env python3
"""Regenerate examples/corsica/corsica.tif from the WTM island_equilibrium clip.

This is the *same* Corsica clip and ocean processing the Water Table Model (WTM)
uses for its island example, so the DepressionHierarchy this repo builds is
directly comparable to the depressions WTM's Fill-Spill-Merge fills with lakes.

Source clip: models/WTM/examples/island_equilibrium/corsica_gebco.tif
  156 x 240 cells, GEBCO 30" (lon 8.4-9.7, lat 41.2-43.2), native NoData = NaN.

WTM's processing (demo.py make_corsica), reproduced verbatim:
  land  = cells where dem > 0, floored to max(dem, 0.5)
  ocean = 0.0 (a flat, unambiguous sea level -- NOT NaN; see ENH-6)
  a 1-cell ocean ring is forced on every edge so the island is fully enclosed.

Run the DH on the result with ocean level 0:
  ./build/dephier.exe examples/corsica/corsica.tif out/corsica 0
"""
import os
import numpy as np
import rasterio
from rasterio.transform import from_bounds

WTM_TIF = os.environ.get(
    "CORSICA_GEBCO",
    os.path.expanduser("~/models/WTM/examples/island_equilibrium/corsica_gebco.tif"),
)
OUT = os.path.join(os.path.dirname(__file__), "corsica.tif")

with rasterio.open(WTM_TIF) as s:
    dem = s.read(1).astype("float32")
H, W = dem.shape
mask = (dem > 0).astype("float32")
mask[0] = mask[-1] = mask[:, 0] = mask[:, -1] = 0          # forced ocean edge ring
topo = np.where(mask > 0, np.maximum(dem, 0.5), 0.0).astype("float32")   # land>=0.5, ocean=0

with rasterio.open(OUT, "w", driver="GTiff", height=H, width=W, count=1,
                   dtype="float32", crs="EPSG:4326",
                   transform=from_bounds(0, 0, W, H, W, H)) as o:
    o.write(topo, 1)

print(f"wrote {OUT}: {W}x{H}, {int((mask > 0).sum())} land / {int((mask == 0).sum())} ocean cells")
