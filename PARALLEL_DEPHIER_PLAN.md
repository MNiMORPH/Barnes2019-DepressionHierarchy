# Implementation Plan: distributed-memory DepressionHierarchy build

**Date:** 2026-07-26
**Status:** PLAN. Companion to `PARALLEL_DEPHIER_DESIGN.md` (the *why*); this is the *how*.
**Repo:** fork `MNiMORPH/Barnes2019-DepressionHierarchy` (upstream `r-barnes/…`).
**Model:** Barnes (2016), *Parallel Priority-Flood* — tile locally, reconcile at boundaries —
generalized from reconciling fill **levels** to stitching the depression **hierarchy**.

This plan is grounded in a close read of `include/dephier/dephier.hpp` (the serial build),
`include/dephier/DisjointDenseIntSet.hpp` (the union-find), and the existing test/driver code.
Line references below point at the code as it stands on the fork's `master`.

---

## 0. Objective and non-goals

**Objective.** Produce a `DepressionHierarchy` for a DEM too large to hold on one node, by
tiling the grid across MPI ranks, so that per-rank memory is **O(N/P) + O(boundary + #depressions)**
and no rank ever materializes a full-grid array. Correctness measured against the serial build
(see §6 for the *right* notion of "correct").

**Non-goals (v1).**
- Distributed Fill-Spill-Merge traversal (that is the Part-3 / FSM repo; it *rides on* the tree
  this build produces — see design note §3). Out of scope here except that our output format must
  be consumable by a distributed FSM.
- Squeezing compute time. Rebuild cadence is decadal–centennial for WTM (design note §4), so this
  is an amortized cost. **Optimize for fit, not speed.** A serial reconciliation phase is
  acceptable in v1 if it fits in memory.

---

## 1. What the serial build actually does (the anchor)

`GetDepressionHierarchy<elev_t, topo>(dem, label, flowdirs)` (`dephier.hpp:263`) runs four phases.
Characterizing each by **memory** and **grid-touch** is what tells us how to distribute it:

| Phase | Code | Touches full grid? | Memory | Notes |
|------|------|--------------------|--------|-------|
| **A. Seed** | `dephier.hpp:305–413` | yes (2 parallel sweeps) | O(N) | find ocean-border + pit cells; sort seeds for determinism |
| **B. Flood + outlet discovery** | `dephier.hpp:451–564` | yes (radix-heap PQ over all cells) | **O(N)** — DEM, `label`, `flowdirs`, PQ, outlet hash | the memory wall. Grows leaf depressions; records the **lowest saddle between each adjacent depression pair** into `outlet_database` |
| **C. Hierarchy assembly** | `dephier.hpp:591–733` | **no** — operates only on the outlet list + `depressions` + `DisjointDenseIntSet` | O(#outlets + #depressions) | sort outlets by elevation, sweep low→high, union-find to build the binary merge tree; special-cases ocean links |
| **D. Volumes** | `dephier.hpp:749–835` | `CalculateMarginalVolumes` yes (per-cell walk up tree); `CalculateTotalVolumes` no (bottom-up over tree) | marginal O(N) read-only; total O(#depressions) | |

**The load-bearing observation for parallelization:** *Phase C never touches the grid.* It reads
`outlets` (one record per adjacent-depression pair) and `depressions` (one record per depression),
and mutates a `DisjointDenseIntSet` keyed on **depression labels, not cells**. Its size is
O(#depressions) ≪ N — the design note's "light tree." So the serial Phase C code is **reusable
almost verbatim** on a *global* outlet set, on one node, if we can produce that global outlet set.

That reframes the whole build: the expensive, must-distribute part is A+B (and the per-cell half of
D); the genuinely new research is only *how boundary crossings enter the global outlet set* so that
Phase C, run once globally, yields the right tree.

---

## 2. Proposed decomposition

```
        ┌──────────── per rank r (tile r) ────────────┐
  A+B   │ run serial A+B on tile r → local labels[r],  │   memory here: O(N/P)
        │ local depressions[r], local internal outlets │
        │ + extract boundary strips (elev,label) for   │
        │   each shared edge                           │
        └──────────────────────┬───────────────────────┘
                               │  gather (small): boundary strips
                               │  + local depression records + internal outlets
                               ▼
        ┌──────── reconcile (rank 0 in v1) ────────────┐
        │ 1. assign global label namespace (prefix sum)│   memory here: O(#dep + boundary)
        │ 2. discover cross-tile outlets from strips    │
        │ 3. GLOBAL Phase C  (reuse serial sort+union-  │
        │    find) over internal ∪ cross-tile outlets   │
        │    → global depression tree                   │
        └──────────────────────┬───────────────────────┘
                               │  broadcast (small): global tree + label remap
                               ▼
        ┌──────────── per rank r ──────────────────────┐
  D     │ remap local labels → global; CalculateMarginal│   memory here: O(N/P)
        │ Volumes on tile r's cells → partial (cell_    │
        │ count,total_elev) per depression; reduce      │
        └──────────────────────┬───────────────────────┘
                               ▼
              CalculateTotalVolumes on rank 0 (O(#dep))
```

**Design choice — centralize Phase C in v1 (recommended).** Gather the O(#depressions) outlet/record
set to rank 0, run the *existing* serial `sort` + `DisjointDenseIntSet` assembly, scatter the tree
back. Rationale:
- Phase C is already the grid-free "light" part; it is the least of the memory problem.
- Reuses tested serial code → smallest new-research surface, most durable (design principle: reuse
  existing patterns rather than invent).
- Amortized/decadal rebuild → a serial reconciliation is fine on wall-clock grounds.
- **Gate:** this assumes the global outlet+depression set fits on one node. That is exactly the
  boundary-graph bound the design note flags as load-bearing (§4). **Measure it first** (§8) before
  committing; distribute Phase C (Barnes 2016 fully-distributed join) only if the bound fails.

This is a cleaner split than "the boundary graph does everything": distribute what is O(N),
centralize what is already O(#depressions).

**Transport.** MPI, one tile per rank, matching Barnes 2016's model; keep the existing OpenMP
parallelism *within* each tile's A/B/D (it is already there — e.g. `dephier.hpp:310, 364, 773`).
Hybrid MPI+OpenMP, as in Barnes 2016.

---

## 3. The research core: cross-tile depression stitching

This is the one genuinely new piece (flood is 2016; DH is 2019; distributed DH is unpublished).

**Boundary strips.** For each shared edge between tiles r and s, each rank contributes the one-cell-
deep strip of `(elevation, local_label)` along that edge (plus the D8 diagonal neighbours at
corners). Size O(tile perimeter) = O(√(N/P)) per edge ≪ N.

**Cross-tile outlet discovery — mirror the serial rule exactly.** In serial, when two differently-
labelled cells are adjacent, the outlet between their depressions is the **higher of the two cells**,
and the database keeps the **lowest such outlet** per depression pair (`dephier.hpp:520–558`). At a
tile boundary we apply the identical rule to each adjacent cross-boundary cell pair `(a∈r, b∈s)`:
candidate outlet elevation `= max(elev(a), elev(b))`; for each depression pair
`(global_label(a), global_label(b))` keep the minimum candidate. This yields cross-tile outlet
records in exactly the `Outlet<elev_t>` form Phase C already consumes. **Merge** these into the
gathered internal-outlet set, then run Phase C once, globally.

### 3.1 The identity question — and why the criterion must be semantic, not structural

The design note (§4) states the distributed tree should be "identical to the serial one." Reading the
flood loop, I believe **exact structural identity does not hold**, and we should not target it. Here is
the reasoning (derived from the code; to be confirmed empirically — see §6):

> Consider a single physical depression (bowl) that straddles the boundary between tiles r and s.
> **Serial:** one wavefront floods it from its single pit cell ⇒ **one leaf** depression, pit at the
> true minimum, outlet at the bowl's rim.
> **Distributed:** each tile floods only its half. Within tile s, the bowl's lowest cells lie *on the
> shared edge*; their lower neighbours are across the boundary (not in s's grid), so they register as
> **pit cells of tile s** (`dephier.hpp:372–388` tests neighbours *in grid* only). So tile s grows its
> own leaf for its half, and tile r grows a leaf for its half. The cross-tile outlet between them sits
> at the col where the bowl crosses the tile line — *below* the rim, *above* both pits. Phase C
> therefore **merges the two halves into a meta-depression** whose outlet is the true rim.

Net: where serial has one node, distributed has three (two boundary-artifact leaves + one meta-node).
**The distributed tree is a refinement of the serial tree**, obtained by inserting internal nodes at
tile crossings. The serial depression still appears — as the meta-node — with the correct rim outlet.

**This does not force a choice between "correct for compute" and "identical for analysis" — §3.2
gives both.** The refinement tree is already correct for every volume/water calculation (below); an
optional post-pass contracts the tile artifacts to recover a serial-identical tree for direct
analysis.

**Consequence for the acceptance criterion (§6):**
1. **Volume conservation (strongest, cheapest check).** For the whole domain and for each serial
   depression, total contained volume matches. `CalculateTotalVolumes` sums children
   (`dephier.hpp:822–825`), so a serial leaf's volume must equal the sum over its distributed
   sub-pieces + the meta-node's marginal. This is a floating-point-tolerance equality and bites hard
   on bugs. Holds on the *refinement* tree, before any collapse.
2. **Outlet-elevation / spill agreement.** Every serial depression's `out_elev` and outlet cell
   matches the corresponding distributed (meta-)node's.
3. **Bit-identity after collapse (§3.2).** After the artifact-contraction pass, the tree is
   node-for-node identical to serial (modulo per-cell `flowdirs` near cuts — see §3.2). This is the
   strong criterion the design note wanted; it is *earned* by the collapse pass, not by the raw build.

### 3.2 Collapsing tile-artifact nodes to a single global DH (Wickert)

The refinement of §3.1 is exactly reversible by a post-processing pass **on the light tree**
(O(#depressions), no full grid — only the boundary strips we already gathered). The build stays
simple and correct-for-compute; this pass, run when a single global DH is wanted *for analysis*,
contracts each split bowl back to one leaf identical to the serial one.

**Collapse criterion — the load-bearing detail.** "The connecting outlet lies on a tile boundary" is
**necessary but not sufficient**: two genuinely distinct bowls whose real saddle happens to fall on the
tile line form a *true* meta-depression that must **not** be collapsed. The distinguishing local
signal comes straight from the pit test (`dephier.hpp:372–388` checks in-grid neighbours only):

> A leaf is a **tile-cut artifact** iff its pit cell sits on a tile edge **and has a strictly-lower
> neighbour across that edge**. A monotonic slope crossing the tile line produces such a spurious pit;
> a real saddle-on-the-line does not (its pits are interior local minima). Collapse only meta-nodes
> whose children reduce to artifact leaves linked by such monotonic crossings; keep meta-nodes over
> real depressions even when their saddle lies on the boundary.

**Defining "strictly-lower neighbour across that edge."** A boundary pit cell lies *on* the tile edge,
so 3 (edge) or 5 (corner) of its 8 D8 neighbours fall in the adjacent tile — the very neighbours the
local pit test skips (`dephier.hpp:378` continues on out-of-grid cells). `cross_neighbour` is one of
those immediately-adjacent cells just across the boundary. Use **strict `<`**
(`dem(cross_neighbour) < dem(pit)`), matching the serial pit definition exactly (`dephier.hpp:380`):
a cell that is a true local minimum stays a pit under any tiling; an artifact pit is one the tile edge
manufactured by hiding a lower neighbour. Strict `<` has **no false positives** — a strictly-lower
cross neighbour means water really would flow across, so the pit is genuinely spurious.

**Deferred residual — flats spanning the boundary.** The strict test does *not* resolve the
**equal-elevation** crossing: a flat-bottomed bowl straddling the line has boundary pits whose cross
neighbours are *equal*, not lower, so neither half is flagged and the single flat stays as two nodes.
This is harmless — volumes still conserve; it is only a spurious extra node in the analysis tree — and
disambiguating it (one flat cut in two vs. two equal-level flats meeting) needs flat-region–aware
analysis, exactly the vast-flats pathology §8 flags. **Ship strict `<` first; treat boundary-spanning
flats as a known, deferred residual** to fold in when flats are handled generally. (`<=` is *not* the
fix — it would over-collapse two genuinely distinct equal-level depressions that merely touch.)

**Algorithm (handles multi-way splits — corners, long boundaries cut by several tiles):**
1. Flag artifact leaves by the pit test above.
2. Form connected components of artifact leaves linked through cross-tile *monotonic-crossing*
   outlets (a corner bowl → 3–4 pieces; a long bowl → a chain).
3. Contract each component and its internal meta-nodes to one leaf `L`:
   `L.pit_cell/pit_elev = argmin/min` over pieces; `L.out_cell/out_elev =` the component's top meta
   outlet (the real rim); `cell_count/total_elevation/dep_vol` are already the whole-bowl totals from
   `CalculateTotalVolumes`, so they carry over unchanged.
4. Rewire references from *outside* the component that pointed at a contracted piece
   (`parent`, `odep`, `geolink`, `ocean_linked[]`) to `L`. O(#depressions) remap.

**Caveat (footnote, not a blocker):** per-cell `flowdirs` near a cut differ from serial (tile-s cells
point to the boundary pit, not across the line). This is irrelevant to pooled-water / FSM behaviour,
but the `flowdirs` field is not serial-identical even after collapse. If a consumer needs
serial-identical flowdirs, that requires a boundary-aware flowdir fix-up — deferred, flagged for
Richard (§10).

This pass is the constructive form of §3.1's contraction criterion: it turns "isomorphic modulo tiling
artifacts" into an executable step and lets the oracle assert **bit-identity** (§6).

---

## 4. Data-structure changes (concrete)

Grounded in the actual structs:

- **Global label namespace.** Each tile's local build makes its own `OCEAN=0` and labels
  `[0, n_r)`. Assign global ids by exclusive prefix sum over `(n_r − 1)` (every tile shares the one
  global OCEAN=0): `global(local=0) = 0`; `global(local=k) = offset_r + (k−1)`. Remap **every label-
  valued field** in each `Depression` record: `parent, odep, geolink, lchild, rchild, dep_label,
  ocean_linked[]` (`dephier.hpp:54–81`) and the boundary-strip labels.
- **Cell indices become global.** `pit_cell` / `out_cell` are `flat_c_idx` (**`uint32_t`**,
  `dephier.hpp:32`) local flat indices. For cross-tile meaning they must map to a global cell
  identity — either `(tile_id, local_flat)` or a global flat index. **Headroom check:** global 30″ is
  N = 933,120,000 < 2³² = 4.29e9, so a global `uint32` flat index *fits* — but only ~4.6× headroom;
  document the ceiling and prefer `(tile,local)` pairs to avoid a silent overflow if resolution grows
  (15″ would blow it). `dh_label_t` is also `uint32` (`dephier.hpp:31`) ⇒ ≤ 2³² depressions globally;
  bound #depressions on real tiles (§8) to confirm headroom.
- **Reuse Phase C unchanged.** `DisjointDenseIntSet` already dynamically grows
  (`DisjointDenseIntSet.hpp:27`) and operates on labels; the global assembly instantiates it at
  `#global_depressions` and runs the loop at `dephier.hpp:640` verbatim. This is the payoff of the
  centralized-C choice.

---

## 5. Public API shape

Keep the serial `GetDepressionHierarchy` untouched (it is the oracle). Add a sibling:

```cpp
template<class elev_t, Topology topo>
DepressionHierarchy<elev_t> GetDepressionHierarchyDistributed(
  TiledArray2D<elev_t>  &dem,      // rank-local tile + halo/edge access + tiling metadata
  TiledArray2D<dh_label_t> &label, // rank-local, filled in place (global labels on return)
  TiledArray2D<int8_t>  &flowdirs, // rank-local
  MPI_Comm comm
);
```

`TiledArray2D` is the new abstraction (tile bounds, neighbour-rank map, edge-strip extraction). It can
wrap richdem's `Array2D` per tile so the inner A/B/D code paths reuse the existing kernels with
minimal edits.

---

## 6. The oracle — which must be *built*, not assumed

**Correction to the design note (§4/§5):** it treats `src/dephier_paper_tests.cpp` as "the repo's
correctness tests / the oracle." In fact that program (built via `CMakeLists.txt:21`) and
`run_tests.sh` form a **crash-only smoke test**: they run the *serial* build on 26 small DEMs
(`test_cases/*.dem`) and write CSV/label outputs, with `set -e` catching only crashes. Nothing
compares outputs to a reference or checks any invariant. **There is no equivalence oracle today.**

**Terminology.** The *oracle* is the **serial `GetDepressionHierarchy` output, treated as ground
truth**. The *differential tester* is the **run that compares serial vs. parallel** on the same DEM:
canonicalize both trees and assert the equivalence criteria below. The two are distinct — one is the
reference, the other the comparison harness — and neither exists in the repo today.

**The test does not need real MPI (the final build does).** Stitching correctness — namespace remap,
cross-tile outlets, collapse — is independent of the transport layer, so the "parallel" side of the
comparison can be an **in-process logical tile-split** (§7.1): one process, tiles split logically, the
exact same stitching code, diffed against serial in plain CI. MPI is required only by the *final*
implementation, where it delivers the actual O(N/P) across-node memory distribution; the MPI build
later re-runs this same comparison over localhost ranks to check the transport didn't break anything.

So a first-class deliverable is the **differential tester**:
1. `serialize(DepressionHierarchy)` → a canonical, order-independent form (sort nodes by
   `(pit_elev, out_elev, cell_count)`; emit outlet elevations and per-node volumes).
2. On each `test_cases/*.dem`: run serial (ground truth) and distributed with an *artificial* tiling
   (e.g. 2×2, 3×1, 2×3, a tiling whose boundary deliberately bisects a known depression, and a
   corner tiling that splits one bowl 4 ways to exercise multi-way collapse).
3. Assert, on the **refinement** tree (pre-collapse): volume conservation to fp tolerance +
   outlet-elevation agreement (§3.1).
4. Assert, on the **collapsed** tree (§3.2): **bit-identity** to the serial tree (node-for-node, all
   fields except per-cell `flowdirs` near cuts). This is the strong check and also validates that the
   collapse criterion neither over-collapses (real saddle-on-boundary) nor under-collapses.
5. Add a degenerate-tiling identity check: a **1×1 tiling must reproduce the serial tree bit-for-bit**
   with no collapse needed — this isolates namespace/remap bugs from genuine stitching effects.

Wire it into CTest: the in-process variant runs everywhere with no MPI at all; the MPI variant runs
the same comparison over a few `localhost` ranks where MPI is available.

---

## 7. Minimal executable tests at each decision (per the numerical-model workflow)

Do *not* build the MPI harness before validating the pieces serially. Order:

1. **Namespace remap, in-process.** Split one small DEM into two logical tiles *in a single process*,
   build each with the serial code, remap into a global namespace, hand the merged outlet set to the
   existing Phase C, and diff against the serial whole-DEM build using the §6 tester. No MPI yet.
   This validates §3 stitching + §4 remap in isolation — a ~100-line harness catching the core
   research risk in seconds.
2. **Boundary bisecting a depression.** Hand-craft a tiny DEM with one bowl astride the split; confirm
   the refinement/volume-conservation behavior predicted in §3.1 (and *discover* it empirically if my
   reasoning is wrong — that discovery is the finding, not a failure).
3. **Then** wrap in MPI (2 ranks, localhost), re-run the tester.
4. **Then** scale ranks and tile counts.

Each step is a runnable check before the next layer of infrastructure — a 10-line test that prevents a
200-line rewrite.

---

## 8. Precursor measurement (do this first — the whole memory argument rests on it)

Before writing distributed code, **bound the reconciliation footprint on real global DEM tiles**
(design note §4, §6 both call this the load-bearing assumption):
- On MERIT/GEBCO-class 30″ tiles, run the *serial* build per tile and count: #depressions per tile,
  #internal outlets, and (by exchanging strips offline) #cross-tile outlets. Extrapolate the global
  #depressions and boundary-outlet totals.
- **Decision gate:** if global (#depressions + #outlets) × record size fits comfortably on one node,
  centralized Phase C (v1) stands. If not, escalate to the distributed join. Report the number
  either way — a "does not fit" result is the finding that redirects the architecture.
- Stress the pathological case the note names: vast flats spanning many tiles (flats become many pit
  cells, `dephier.hpp:356–361`), which could inflate boundary outlets.

---

## 9. Sequencing (revised from design note §5)

0. **Precursor:** bound reconciliation size on real tiles (§8). *Gates the v1 architecture.*
1. **Differential oracle** (§6) — build it against the *serial* code first; it has value immediately
   and is prerequisite for trusting everything after.
2. **In-process tile-split + namespace remap + reuse Phase C** (§7.1) — validates the research core
   with no MPI.
3. **Cross-tile outlet discovery + boundary-bisection semantics** (§7.2, §3.1).
4. **Artifact-collapse pass** (§3.2) — pure tree post-processing; unlocks the bit-identity oracle
   check (§6.4). Can be built and tested against the in-process split (step 2) before any MPI.
5. **MPI harness:** `TiledArray2D`, edge-strip exchange, gather/scatter around centralized Phase C
   (§2, §5).
6. **Distributed Phase D volumes** (§2) with reduction; verify volume conservation at scale.
7. **Scale + load-balance study** (design note §4); refine tiling only if measured imbalance warrants.
8. *(Later / separate repo)* distributed FSM traversal on the produced tree.

---

## 10. Open questions for Richard (co-scoping)

1. **Two trees: refinement for compute, collapsed for analysis (§3.1–§3.2).** Plan is: the build
   yields a *refinement* tree (correct for all volume/water calcs), and an optional O(#depressions)
   post-pass contracts tile-cut artifacts to a serial-identical tree for direct analysis — so we get
   both without a boundary-aware re-flood. Does the collapse **criterion** (artifact leaf = boundary
   pit with a lower cross-edge neighbour; keep real saddles-on-the-line) hold in cases we're not
   thinking of? And do any consumers need serial-identical per-cell `flowdirs` near cuts (the one
   thing collapse does *not* restore)?
2. **Centralized vs. distributed Phase C (§2).** Given decadal rebuild cadence, is a serial
   reconciliation acceptable for v1, deferring the fully-distributed 2016 join until §8 says we must?
3. **Reuse of the 2016 join code.** Is there existing Barnes-2016 join/tile-exchange code (in richdem
   or elsewhere) we should build `TiledArray2D` and the strip exchange on, rather than writing fresh?
4. **Global cell-index convention (§4).** `(tile,local)` pair vs. global `uint32` flat index — any
   downstream (FSM) constraint that forces one?
```
