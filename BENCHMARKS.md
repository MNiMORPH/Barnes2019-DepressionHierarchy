# Benchmarks

## Footprint sweep — distributed flat-label replay (`DH_FLAT_PARTITION_REPLAY`)

**What & why.** The distributed DepressionHierarchy build's design goal is a *footprint-bounded* per-rank
memory profile — `O(N/P)` owned cells plus an `O(boundary)` halo — so it scales to grids no single rank
can hold. The one stage whose footprint is not obviously bounded is the flat-label replay
(`DH_FLAT_PARTITION_REPLAY`), which reproduces serial's exact flat partition across seams by growing an
**adaptive halo** until each owned cell's replay label is stable. A cell's replay label only settles once
the halo reaches that cell's basin *pit* — so on real terrain the uncapped halo grows to *basin* extent,
which can be most of the grid. A finite **halo cap** (the 4th CLI argument) bounds it to
`owned_cols + 2·cap`, at the risk of no longer reproducing serial's exact tree inside basins wider than
the cap (the result is still a valid flow field and tree). This sweep quantifies both axes.

**Metric.** With `DH_HALO_DIAG=1`, each rank prints `owned_cols` (`= W/ntiles`, the `O(N/P)` baseline) and
`held_cols` (owned + halo = the peak transient footprint of the replay). We report the **peak `held_cols`
across ranks**, plus whether the tree stayed bit-identical to serial (`MPI-TREE-MATCH`) and the flowdir
diff. Wall-clock is deliberately *not* measured: the thread shim is not a performance vehicle (serialized
message queue, threads ≠ processes), and the fixtures are tiny — footprint is the meaningful axis here.

**Reproduce.** `tools/footprint_sweep.sh` (needs the thread-shim `build/dephier_mpi.exe`; by default
generates a 256×256 β=1.5 fractal via `tools/make_synthetic_dem.py`, which needs `WTM_DIR` — or pass your
own DEM: `tools/footprint_sweep.sh my.dem`).

### Results — 256×256 fractal (β=1.5, seed 1), `W=256`

**Sweep A — halo cap at fixed `ntiles=8`** (`owned_cols=32`):

| cap | held_cols (peak) | held/W | tree | flowdir |
|-----|------------------|--------|------|---------|
| ∞ (default) | 224 | 0.88 | MATCH | fd_diff=0 |
| 64 | 160 | 0.62 | MATCH | fd_diff=0 |
| 32 | 96 | 0.38 | MATCH | fd_diff=0 |
| 16 | 64 | 0.25 | MATCH | fd_diff=0 |
| 8  | 48 | 0.19 | MATCH | fd_diff=0 |
| 4  | 40 | 0.16 | **DIFFER** | fd_diff=0 |

**Sweep B — `ntiles` at fixed `cap=16`** (footprint = `owned_cols + 2·cap`):

| ntiles | owned_cols | held_cols (peak) | tree | flowdir |
|--------|-----------|------------------|------|---------|
| 2 | 128 | 144 | MATCH | fd_diff=0 |
| 4 | 64  | 96  | MATCH | fd_diff=0 |
| 8 | 32  | 64  | MATCH | fd_diff=0 |

### Takeaways

- **Uncapped replay is exact but not footprint-bounded.** The default (∞) halo grew to 224/256 columns
  (88% of the grid) even on ordinary fractal terrain — a single rank ends up holding most of the grid.
- **A finite cap bounds the footprint linearly:** `held_cols = owned_cols + 2·cap`, confirmed exactly by
  Sweep B (e.g. 32 + 2·16 = 64). Per-rank footprint is `O(N/P) + O(cap)`.
- **Exactness has a per-DEM cap threshold.** On this DEM the tree stayed **bit-identical** down to `cap=8`
  (holding only 19% of `W`) and first diverged at `cap=4`. So a modest cap buys a ~4–5× footprint cut at
  *zero* exactness cost — a cap of ~8–16 is a good default here, and worth setting rather than leaving
  the halo unbounded.
- Flowdirs are unaffected by the flat-replay cap (`fd_diff=0` throughout — flowdirs come from the separate
  flat-flowdir pass).

The right cap is DEM-dependent (it scales with the widest basin that must be reconciled across a seam);
the sweep is the way to pick it for a given tiling and terrain.
