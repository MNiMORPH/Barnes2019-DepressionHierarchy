# Notes for Richard — changes to the DepressionHierarchy core

**Context:** distributed-memory DH build on the MNiMORPH fork (see `PARALLEL_DEPHIER_DESIGN.md`,
`PARALLEL_DEPHIER_PLAN.md`). This file collects the changes to *your* core code
(`include/dephier/dephier.hpp`, `include/dephier/radix_heap.hpp`) that warrant your review before any
upstream merge — especially the ones that **change serial output**. Nothing here has been pushed
upstream; it lives on the fork's `master`.

**PR status (2026-07-30):** Andy has raised these changes with Richard, who is amenable to them going
upstream — this file is the review note for that PR. Two things to sharpen before opening it: (a) confirm
which of the serial-output-changing items below (#1 `PhaseCD` tiebreak, #2 `radix_heap` order, #3 `out_cell`
tie-break) Richard has explicitly OK'd vs. is still reviewing; (b) note the intended perf follow-up on #2
(the composite radix key, so the flat sort is not a per-bucket cost). Only the `dephier.hpp` /
`radix_heap.hpp` changes (the two lists below) are candidates for the upstream PR; the distributed harness
(`dephier_mpi.cpp`, `dephier_stitch.cpp`) is fork-only. Status since the last note: the distributed build is
complete (ENH-8) and **bit-identical to serial** across the fixture sweep — that bit-identity rests on the
serial-determinism changes #1–#3 here **plus** an exact *reproduction* of your flat labelling (the pit-index
flood replay — NOT a serial change; it matches your algorithm rather than altering it).

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

3. **`PhaseAB` outlet `out_cell` tie-break** (`dephier.hpp`, outlet discovery — the last split-invariance
   piece; 2026-07-29). When two connections between the same depression pair share the lowest `out_elev`,
   the kept `out_cell` was whichever the flood reached FIRST (pop order) — traversal-order-dependent, so a
   tiled build kept a *different* `out_cell` than the whole-grid flood on nearly every tied outlet. Since
   your `PhaseC` sort breaks `out_elev` ties by `out_cell` (#1), that flips the tree at tied outlets → an
   equally-valid but different meta hierarchy. Fix (one line): on an `out_elev` tie keep the **lowest-index**
   `out_cell` — geometric, traversal-order-independent, the same rule the stitch's re-derived outlet scan
   already uses. Your `Depression::out_cell` comment already called this tie "arbitrarily chosen"; this makes
   the choice canonical, and makes **serial itself tiling-independent at tied outlets**.
   **Effect (measured):** with this change *and* our exact reproduction of your flat labelling (see the
   residual section — no serial change needed there), the tiled build is **bit-identical to serial on all
   107 fixture×split cases** (was 79/107). Argument, not just data points: once both sides feed `PhaseC` the
   same outlets and the leaf partition already matches, the deterministic sort yields identical trees by
   construction. Our regression suite (assertions are serial-relative) still passes 25/25.

## Changes that PRESERVE serial output (FYI, verified byte-identical)

4. **`GetDepressionHierarchy` split into exposed phase functions** (`dephier.hpp`).
   `FloodAndAssignDepressions` = flood + outlet discovery (paper A-C), exposing `{depressions, outlets}`;
   `ConstructHierarchy` = the grid-free hierarchy assembly (paper D); `ConstructHierarchyAndVolumes` =
   `ConstructHierarchy` + volumes (§6.4), the central-only convenience. `GetDepressionHierarchy` is retained
   as a thin serial wrapper (unchanged signature). A distributed build runs `FloodAndAssignDepressions` per
   tile, one global `ConstructHierarchy`, and distributes the volume step. (Named to the paper, not shifted
   letters — see `DH_API_NAMING_REVIEW.md`; earlier entries above still cite the old `PhaseAB/PhaseC/PhaseCD`
   names, mapped there.)
   Verified: canonical signature identical before/after across all `test_cases` + synthetic trees, with
   a negative control (a deliberate perturbation moved 26/28 signatures).

5. **`BOUNDARY` exterior label** (`dephier.hpp`).
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

   **→ RESOLVED (2026-07-29), and NOT by changing your labels — we reproduce them exactly.** Rather than a
   new flat-LABEL rule, a full pit-up replay of *your* flood's label partition (seeds = ocean + every
   cell with no strictly-lower neighbour; the radix `(elevation, index)` pop order) is provably
   partition-identical to serial on 49 DEMs (all fixtures + Corsica + adversarial fractals). So the
   cross-seam flat-label divergence closes by *matching* your algorithm, not altering it — no serial-label
   change. (Today a full-grid pass behind a flag in the stitch; the distributed 1-column seam-exchange form,
   replaying flood order across the seam, is engineering on our side, still not a serial change.)

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

   **→ RESOLVED (2026-07-29) by change #3 in the serial-changing list above.** The diagnosis here is exactly
   right — "serial's flood breaks the outlet-cell tie by flood order, the geometric re-derivation by lowest
   index, and those differ." The resolution is to make *serial* break it by lowest index too (the one-line
   `PhaseAB` change), so both sides agree without the re-derivation having to reproduce flood order. With that
   plus the flat-label reproduction (#2), the tiled build is bit-identical to serial on all 107 cases. This
   *was* the whole residual: every remaining DIFFER, audited, had matching pairs + `out_elev` and differed
   only in `out_cell`.

## Also for your review (shared-core input handling — ENH-6)

`ocean_labels` treats every NoData cell as OCEAN, and a land↔ocean outlet's sill is taken from the ocean
cell's **raw DEM value** (the NoData sentinel). So the tree depends on the arbitrary sentinel: `9` "works"
by luck, `-9999` shifts sills, **`NaN` aborts** `PhaseC`'s outlet sort (real GEBCO NoData). An OCEAN cell
is base level — the sill should be the **land** cell's elevation, and NoData's meaning (ocean vs interior
void) should be declared, not guessed. Proposed: a base-level sill in the shared outlet code + a
`--nodata {ocean|void|error}` flag. Touches `ocean_labels` / outlet elevation, hence your bucket. (Tracked
as ENH-6, GitHub issue #2.)
