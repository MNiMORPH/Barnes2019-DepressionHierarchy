# Benchmarks

## Choosing a mode: footprint vs. exactness

The distributed build has three run modes, trading per-rank memory footprint against how closely the
result reproduces the serial tree. **Pick by whether the grid fits in one rank's memory and whether you
need bit-identity — not by defaulting to "most exact".** The flat-label replay is *opt-in*
(`DH_FLAT_PARTITION_REPLAY`); the halo cap is the optional 4th CLI argument.

| mode | invoke | tree vs. serial | per-rank footprint | use when |
|------|--------|-----------------|--------------------|----------|
| **default** (no replay) | (flag off) | **volume-correct + valid**, not bit-identical (flat cells may fall in a different valid leaf) | **bounded**: `O(N/P) + O(boundary)` (1-column halo) | **production at scale** — grids too big for one rank; volume-correctness is enough |
| **uncapped replay** | `DH_FLAT_PARTITION_REPLAY=1`, no cap | **bit-identical** | **unbounded**: the adaptive halo grows to basin extent, ~O(W) — a rank can hold most of the grid | you need bit-identity **and** the grid fits in one rank (validation; small/medium DEMs) |
| **capped replay** | `…=1` + cap arg | valid + volume-correct; bit-identity is *jagged* in the cap (see below) | **bounded**: `held_cols = owned_cols + 2·cap` | you want maximal exactness **and** it must fit at scale — the narrow niche the cap exists for |

**The trap to avoid: "just always use uncapped."** Uncapped is exact, but it re-introduces the O(N)
per-rank footprint the distributed design set out to avoid — on a large, flat-heavy grid (e.g. GEBCO 30″
with the Caspian / Great Lakes) a rank's halo can swell toward the whole grid and **OOM the rank**,
defeating the point of distributing. Uncapped is the right choice *only when memory allows*.

**Rule of thumb:**
- Grid fits in one rank, or you need bit-identity → **uncapped replay** (or just run serial).
- Grid does not fit and volume-correctness suffices → **default (no replay)** — the scalable path.
- Grid does not fit and you want near-exactness → **capped replay**, accepting that a finite cap gives a
  valid, volume-correct tree whose *exact* signature (and, under an aggressive cap, node count) is not
  guaranteed. Only `VOL-MATCH` is guaranteed; see the sweep below for how to pick a cap for a given DEM.

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

### Sweep C — exactness vs. cap across terrain roughness (β), `ntiles=8`

Whether a *capped* run reproduces serial's tree bit-for-bit (`MPI-TREE-MATCH`) vs. produces a valid but
different tree (`MPI-TREE-DIFFER`), swept over spectral roughness β (higher β = smoother/larger-scale
features). `held_cols` is identical to Sweep A (terrain-independent: `owned + 2·cap`). **Every `DIFFER` cell
below is still `MPI-VOL-MATCH`** — a valid, volume-correct tree, never a wrong or degraded one; that is the
*universal* capped-mode invariant. The `DIFFER` cells in *this* sweep are further `MPI-DECOMP-CORRECT` (same
depression count, so the difference is only the fine tie-break signature) — but that stronger property is
**not** guaranteed in general: a sufficiently aggressive cap can change the node count too (see the takeaway).

| β (256², seed 1) | ∞ | 64 | 32 | 16 | 8 | 4 | 2 |
|------------------|---|----|----|----|---|---|---|
| 1.0 (rough)      | M | M | M | M | M | M | **D** |
| 1.5              | M | M | M | M | M | **D** | D |
| 2.0              | M | M | M | M | **D** | M | M |
| 2.5 (smooth)     | M | M | M | M | M | M | M |

Read the β=2.0 row carefully: `cap=8` **DIFFER** but `cap=4` **MATCH** — deterministic (5/5 each).

### Takeaways

- **Uncapped replay is exact but not footprint-bounded.** The default (∞) halo grew to 224/256 columns
  (88% of the grid) even on ordinary fractal terrain — a single rank ends up holding most of the grid.
- **A finite cap bounds the footprint linearly and terrain-independently:** `held_cols = owned_cols + 2·cap`,
  confirmed exactly by Sweeps B and C. Per-rank footprint is `O(N/P) + O(cap)`.
- **Exactness is NON-MONOTONIC in the cap, and has no terrain-independent safe threshold.** A *larger* halo
  is not reliably "more serial-identical": β=2.0 gives the exact serial tree at `cap=4` but a different
  (valid) one at `cap=8`. Different caps resolve different seam-straddling basins, so bit-identity is a
  jagged function of the cap × terrain × tiling. (A single-terrain sweep can look like a clean threshold —
  Sweep A does — but that is a coincidence of that terrain, not a rule.) **It is not a correctness risk:**
  every capped result stays `VOL-MATCH` (volume-correct, valid tree), so the non-monotonicity is movement
  *among valid trees*, never between a correct tree and a broken one. (In this sweep those valid trees are
  also `DECOMP-CORRECT`, so only the tie-break signature moves; a more aggressive cap can additionally change
  the node count — still `VOL-MATCH` — landing in the documented meta-vs-`ocean_linked` tie-break class, e.g.
  kerry_test10 @ 5 tiles cap=1 is `DECOMP-INCORRECT` + `VOL-MATCH`.) The mechanism: the capped replay makes a
  *windowed* geodesic partition of a seam-straddling flat, and the radix-pop-order tie-break at the window
  *edge* is discontinuous in where the edge falls; only a window containing the whole straddling basin
  reproduces serial's exact pop order (which is exactly why the full halo is needed for bit-identity).
- **The cap is a footprint/validity knob, not a footprint/exactness dial.** For guaranteed bit-identity,
  use *uncapped* (or, if you happen to know it, a cap ≥ the widest seam-straddling basin). A finite cap
  always yields a **valid, volume-correct tree** (`VOL-MATCH`) at a bounded footprint — usually with serial's
  depression count too, though an aggressive cap can change that — just do not expect its exact signature (or
  even its node count) to converge monotonically to
  serial as you raise the cap.
- Flowdirs are unaffected by the flat-replay cap (`fd_diff=0` throughout — flowdirs come from the separate
  flat-flowdir pass).
