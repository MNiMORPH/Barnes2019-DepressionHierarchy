# Design Note: a distributed-memory (parallel) DepressionHierarchy build

**Date:** 2026-07-26
**Status:** DESIGN / strategic scoping. For deciding direction, not a commitment.
**Authors:** Andy Wickert + Claude (for discussion with Richard Barnes)
**Context:** motivated by the Water Table Model (WTM; Callaghan & Wickert), whose path to
**global 30″** runs is now blocked on memory, not compute — and the DepressionHierarchy
*build* is the binding full-grid structure. See the companion WTM note
`GLOBAL_SCALING_DESIGN.md` in the WTM repository for the downstream framing.

## 1. Why: it's a memory problem

WTM's groundwater solve is now distributed and cheap, so global WTM is feasible on *compute*.
What it can't do is **fit**: at 30″ the domain is `N = 933,120,000` cells (43200 × 21600), and
a single `double` full-grid field is ~7.5 GB. WTM has distributed its own solve arrays; the
**last, largest full-grid structure is the DepressionHierarchy build**, which this library
performs serially on one node.

The key distinction (and why this note is about the *build*, not the traversal):

- **`GetDepressionHierarchy` (`dephier.hpp`) is memory-heavy** — it holds *several* O(N)
  full-grid arrays: the DEM, the per-cell depression **labels**, flowdirs, plus the
  priority-flood frontier (`radix_heap`) and the `DisjointDenseIntSet` union-find over cells.
- **The resulting hierarchy (the tree) is light** — O(#depressions) ≪ N — and downstream
  traversal (Fill-Spill-Merge, the Part-3 repo) works on that tree plus the water field, which
  a coupled model already holds distributed.

So for global-scale use, **the build is the ceiling**; distributing it is the enabling step.
Reason = memory; any compute speedup from parallelism is a welcome bonus, not the motivation.

## 2. The build already *is* a priority-flood + a merge tree

`GetDepressionHierarchy` is a priority-flood (`priority_flood.hpp`, `radix_heap.hpp`) that, as
the flood rises from the ocean/edges, **records merge events**: when two sub-depressions meet
at a saddle they are unioned (`DisjointDenseIntSet`) and a parent depression is created,
assembling the tree bottom-up. That is precisely the structure that makes a distributed version
tractable:

> **distributed DH build = distributed priority-flood, generalized to carry the merge tree.**

And distributed priority-flood is a *solved problem*: **Barnes (2016), "Parallel Priority-Flood
depression filling for trillion cell digital elevation models on desktops or clusters"**
(*Computers & Geosciences*). This note is that method, extended from reconciling fill *levels*
to stitching the depression *hierarchy*.

## 3. The shape

- **Tile the DEM across ranks** (one subdomain per rank). Each rank runs the existing
  `GetDepressionHierarchy` on its tile, producing its **local** labels, local union-find, and
  local depression tree. *This is where the O(N) memory distributes:* DEM + labels → O(N/P)
  per rank. This step is close to the current serial code, called per tile.
- **Stitch across tiles via a boundary graph.** Build a global structure over only the
  tile-edge spill points (Barnes 2016's "join"). It must do two things:
  1. reconcile cross-tile fill levels / spill directions (the 2016 flood already does this), and
  2. **unify tile-spanning depressions and record boundary saddle/merge events into one global
     hierarchy** — i.e. a `DisjointDenseIntSet`-style union-find promoted to operate over
     *boundary depression IDs* rather than cells. This tree-stitch is the **new** piece (a
     distributed DH has not been published; the flood is 2016, the DH is 2019).
  The boundary graph is size ~O(boundary) **≪ N** — the only whole-domain object any rank holds.
- **Downstream traversal (FSM)** then runs on the distributed hierarchy: tile-local water
  movement with cross-tile spill/merge exchanged through the same boundary graph.

**Memory outcome:** per-rank footprint → **O(N/P) + O(boundary)** → no rank holds a full-grid
array → global-scale DH becomes feasible.

## 4. Risks / open questions

- **Boundary-graph size — the load-bearing assumption.** Everything rests on the inter-tile
  structure being ≪ N. Barnes 2016 already floods *trillion-cell* DEMs, strong evidence it is;
  the DH stitch adds only sparse tree links (~O(#cross-tile depressions)). Still worth bounding
  against real global DEM tiles for pathological cases (vast flats spanning many tiles).
- **Cross-tile hierarchy stitching — the research core.** Correctly promoting the per-cell
  `DisjointDenseIntSet` merges to a *global* union-find over boundary depressions, and recording
  the boundary saddles as parent nodes, so the distributed tree is identical to the serial one.
  Likely needs iteration to global consistency; bound the number of rounds vs P.
- **Determinism / correctness.** The distributed build should reproduce the serial hierarchy
  (the repo's correctness tests, `dephier_paper_tests.cpp`, are the oracle) — a strong,
  testable acceptance criterion.
- **Rebuild cadence (for dynamic-topography users).** WTM rebuilds the DH only
  decadal–centennial (isostasy), so the distributed build is an amortized, occasional cost —
  optimize for *fit*, not *speed*.
- **Load balance.** Depressions cluster spatially; a uniform tiling may imbalance the *build*.
  Barnes 2016's tiling already targets this; keep it simple first, measure, then refine.

## 5. Sequencing

1. **Wrap the existing serial build to run per-tile** (local labels + local hierarchy) — small.
2. **Boundary graph: fill-level/spill reconciliation** (port Barnes 2016 join) — known.
3. **Boundary graph: hierarchy stitch** (global union-find over boundary depressions) — the new
   research; validate against `dephier_paper_tests.cpp` on tiled-vs-serial equivalence.
4. **Distributed FSM traversal** on the hierarchy (in the FSM/Part-3 repo) — rides on the above.

## 6. Bottom line

The DepressionHierarchy *build* — not the traversal — is the memory wall for global-scale
users, and the build is a priority-flood plus a merge tree, so a distributed build is
**Barnes (2016) generalized to stitch the hierarchy across tiles**. It is justified by memory
alone (compute a bonus), it is a research-grade effort, and it is naturally **Richard Barnes'
territory as author of both this library and the parallel priority-flood** — so the intent is
to co-scope and co-author, starting from a fork (MNiMORPH). Recommended precursor: bound the
boundary-graph size on real global DEM tiles, since the whole memory argument rests on it —
concrete input on hand is the **GEBCO 30″ global DEM** (`gebco_08.nc`, the full 933,120,000-cell
grid): tile it and count the spill points / cross-tile depressions against `N`.
