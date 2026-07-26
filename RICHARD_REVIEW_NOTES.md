# Notes for Richard — changes to the DepressionHierarchy core

**Context:** distributed-memory DH build on the MNiMORPH fork (see `PARALLEL_DEPHIER_DESIGN.md`,
`PARALLEL_DEPHIER_PLAN.md`). This file collects the changes to *your* core code
(`include/dephier/dephier.hpp`, `include/dephier/radix_heap.hpp`) that warrant your review before any
upstream merge — especially the ones that **change serial output**. Nothing here has been pushed
upstream; it lives on the fork's `master`.

## Changes that CHANGE serial output (please review)

Both are arguably fixes for *latent non-determinism* — the previous output depended on hash-map /
insertion order, so it could already differ across platforms/STL versions. But they do change the
tree/labels for the affected cases, so they're yours to bless.

1. **`PhaseCD` geometric tiebreak** (`dephier.hpp`, outlet sort).
   The outlet sort broke `out_elev` ties by discovery order (unordered_map iteration). Now it breaks
   ties by **outlet cell, then the endpoint depressions' pit cells** — purely geometric, label-
   namespace-independent. This makes the binary hierarchy deterministic (needed so a distributed build
   reproduces the serial tree at triple junctions, where several outlets share `out_elev` *and*
   `out_cell`). **Effect:** re-binarizes simultaneous merges; ~12 of 28 test DEMs change to an
   equally-valid tree.

2. **`radix_heap` deterministic flat order** (`radix_heap.hpp`, `pair_radix_heap::pull`).
   Equal-elevation items were popped LIFO (insertion order), so **flats resolved differently between a
   whole-grid and a tiled flood**. Now each minimum-key bucket is sorted by value (cell index) so
   equal-elevation cells pop in a fixed geometric order. Row-major index preserves order between global
   and tile-local indexing within a tile, so a within-tile flat resolves identically serial vs tiled.
   **Effect:** changes serial flat resolution to a deterministic one.
   **Perf caveat:** this sorts every equal-elevation bucket — a cost on flat-heavy DEMs. The intended
   optimization is a **composite radix key** (fold the cell index into the key so the radix orders by
   `(elevation, index)` with no per-bucket sort). Not yet done.

## Changes that PRESERVE serial output (FYI, verified byte-identical)

3. **`GetDepressionHierarchy` split into `PhaseAB` + `PhaseCD`** (`dephier.hpp`).
   `PhaseAB` = flood + outlet discovery, exposing `{depressions, outlets}`; `PhaseCD` = hierarchy
   assembly + volumes. `GetDepressionHierarchy` is retained as a thin serial wrapper (AB then CD) with
   the original signature. A distributed build runs `PhaseAB` per tile and one global `PhaseCD`.
   Verified: canonical signature identical before/after across all `test_cases` + synthetic trees, with
   a negative control (a deliberate perturbation moved 26/28 signatures).

4. **`BOUNDARY` exterior label** (`dephier.hpp`).
   New sentinel `BOUNDARY = max-1` for cells that drain off a tile edge in a distributed build. It
   seeds/spreads through the flood exactly like the ocean (so the flood loop is unchanged) but is NOT a
   terminal sink — `PhaseCD` never sees it (the stitch reconciles it across tiles first). `PhaseAB`'s
   three exterior checks were generalized `OCEAN` → `OCEAN||BOUNDARY` (seed gathering, no-seeds guard,
   pit-finding skip); input validation now permits `NO_DEP`/`OCEAN`/`BOUNDARY`. A serial run has no
   `BOUNDARY` cells, so serial output is unchanged (verified byte-identical).

## Open design questions (co-scoping — full list in PARALLEL_DEPHIER_PLAN.md §10)

- **Refinement vs. collapse:** the distributed tree is a *refinement* of the serial one around tile
  boundaries; an O(#depressions) collapse pass recovers a serial-identical tree for analysis. Agree
  that's the right model (vs. targeting bit-identity in the raw build)?
- **Centralized vs. distributed `PhaseCD`:** `PhaseCD` is grid-free (O(#depressions)); v1 centralizes
  it. Acceptable given decadal rebuild cadence, deferring the fully-distributed 2016 join until the
  §8 footprint measurement says otherwise?
- **Reuse:** we intend to build the stitch on your `A2Array2D`/`Layoutfile` tiling + `HandleEdge` join
  (`submodules/richdem/programs/parallel_priority_flood/main.cpp`). Right base, or a newer layer?
- **Seam-straddling flats / halo:** your perimeter-strip join carries no halo. For the *hierarchy*, a
  flat/tie or basin straddling a seam needs cross-tile reconciliation. Any preference from the filling
  work?
