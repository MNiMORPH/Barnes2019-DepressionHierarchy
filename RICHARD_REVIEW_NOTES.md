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

## Remaining split-invariance residual — characterization + one question (2026-07-29)

We raised the bar from "correct-volume + valid-tree" to *split-invariance* (identical tree regardless of
how the array is tiled) and chased the residual. **Every remaining DIFFER is volume-correct and a valid
tree** (`VOL-MATCH` across the sweep) — so under the refinement/collapse model above, none is a bug. But
characterizing them showed the old "PhaseCD tie-break / out_cell class" was **not one thing**; it splits
three ways, and only one is really yours:

1. **A stitch re-derivation bug — FIXED on our side, was never yours.** The stitch re-derives outlets from
   the resolved label grid and lacked the explicit lower-`out_cell` tie-break your `PhaseCD` path already
   had, so its scan was order-dependent (intra-pairs then seam-pairs). Fixed (commit `2948000`): sweep
   MATCH 72→75, `kerry_test4` splits 3/8/9 DIFFER→MATCH. We had *mis-filed* these as "yours." (Consolidated
   all three outlet scans into one shared helper so this can't drift again — ENH-5, `b110f0d`.)

2. **Cross-seam flat LABEL assignment (the one that may be yours).** On nested-flat fixtures the tiled
   flood assigns flat *cells* to depressions differently than the whole-grid flood — `kerry_test2` split 3
   has 12 raw-label diffs, `kerry_test11` split 7 has 407 — which re-binarizes meta-vs-`ocean_linked` at
   tied outlets (same depressions, same volume). Your radix fix (#2 above) made flat *pop-order*
   deterministic so *within-tile* flats resolve identically; but *which depression* a flat cell is labelled
   into across a seam is still flood-order-dependent. **Question:** do you want a geometry-deterministic
   flat-LABEL rule (the label analog of #2 / of the resolve_flats flat-routing you already bless for
   flowdirs), which would make flat labelling split-invariant but **change serial labels** on flat-heavy
   DEMs? Or is the collapse/refinement model sufficient here and we leave it?

3. **Outlet-order meta-vs-`ocean_linked` with MATCHING labels — INVESTIGATED; it folds into #2, not #1.**
   A few cases (`kerry_test3`/`kerry_test7` split 3) have **0 raw-label diffs** — labels bit-identical to
   serial — yet `PhaseCD` builds a meta where serial builds an `ocean_linked` (a big elevation tie: three
   flat stripes whose separating walls and the ocean frame are all at the same elevation). We hoped this was
   a fixable re-derivation tie like #1. It is not. Dumping both outlet sets (`kerry_test7` split 3): serial
   has 10 outlets, the re-derivation 5, and they sit on **different cells** — serial's elevation-6 saddles
   land at high rows (`(3,6)`,`(3,7)`,`(1,7)`), the re-derivation's at low rows (`(3,2)`,`(5,2)`,`(1,1)`).
   The whole rim/wall is a flat at elevation 6, so *which* cell is "the" outlet between a pair is a tie that
   serial's flood breaks by flood order and the geometric re-derivation breaks by lowest index — and those
   differ. So this is the **same cross-seam flat / flood-order root as #2**, manifesting in outlet-CELL
   selection rather than label assignment; the geometric re-derivation cannot reproduce serial's choice
   without reproducing flood order. It is subsumed by the #2 question (a geometry-deterministic flat rule),
   **not** separately fixable on our side. (Correcting an earlier, too-optimistic note that called it fixable
   like #1 — the ENH-5 lesson cutting the other way.)

## Also for your review (shared-core input handling — ENH-6)

`ocean_labels` treats every NoData cell as OCEAN, and a land↔ocean outlet's sill is taken from the ocean
cell's **raw DEM value** (the NoData sentinel). So the tree depends on the arbitrary sentinel: `9` "works"
by luck, `-9999` shifts sills, **`NaN` aborts** `PhaseC`'s outlet sort (real GEBCO NoData). An OCEAN cell
is base level — the sill should be the **land** cell's elevation, and NoData's meaning (ocean vs interior
void) should be declared, not guessed. Proposed: a base-level sill in the shared outlet code + a
`--nodata {ocean|void|error}` flag. Touches `ocean_labels` / outlet elevation, hence your bucket. (Tracked
as ENH-6, GitHub issue #2.)
