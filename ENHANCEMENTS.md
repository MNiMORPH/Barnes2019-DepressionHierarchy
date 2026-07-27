# Enhancements (deferred, designed) — distributed DepressionHierarchy

Local tracker for designed-but-deferred work on the parallel/distributed DH build, so the
plan stays on a linear path to completion. Each entry is written like a GitHub issue and can
be posted to the MNiMORPH fork's tracker when wanted.

---

## ENH-1: bit-identical O(boundary) distributed flat resolution ("option 2")

**Status:** designed, deferred (the problem is already *bounded* by option 3; this recovers
bit-identity on the giant-flat tail without growing per-rank memory). Build on the MPI harness.

**Type:** enhancement / performance-quality. Not a blocker.

### Context
`resolve_flat_flowdirs_distributed` (`tools/dephier_stitch.cpp`) resolves flat flowdirs per tile
with an **adaptive halo** (PARALLEL_DEPHIER_ENGINEERING.md §6, move 2), capped by **option 3** so
footprint is `O(N/P)+O(cap·boundary)` for any input. Consequences of the cap: a flat wider than the
cap is resolved with a capped halo — a *valid* convergent flow field (same tree, same sinks) but not
necessarily serial-identical in that flat's deep interior. Measured flat extents on GEBCO 30″
(`dephier_flat_extent`) top out ~851 cells (Great Lakes), so with large MPI tiles the cap rarely or
never bites; this enhancement matters only when tiles are smaller than the largest lakes, or when
strict bit-identity of the `flowdirs` field is required there.

### What it does
Replace the *capped* branch (only) with the strict-`O(boundary)` **iterative 1-cell gradient
exchange**: reproduce Barnes-2014 `resolve_flats` distributed, so the flat routing is **bit-identical
to serial** with per-rank memory `O(N/P)+O(boundary)` and cost paid in communication *rounds*
(≈ flat diameter in cells), not memory. This is the memory-for-latency trade that suits a
memory-bound problem.

### Design (from reading `richdem/flats/flat_resolution.hpp`)
Barnes-2014's mask is `flat_mask(c) = flat_height[label] − away_dist(c) + 2·towards_dist(c)`, where
`away_dist` = BFS distance from a flat's **high edges**, `towards_dist` = BFS distance from its **low
edges**, and `flat_height[label]` = the **max** `away_dist` in that flat. Distributing it bit-identically
needs four cross-tile pieces — the same machinery the MPI build already provides:

1. **Flat-label unification across seams** — connected flats spanning tiles share one global label
   (union-find over seam-adjacent equal-elevation flat cells; mirrors the depression-label remap, §4).
2. **Two distributed multi-source BFSs** (`away`, `towards`) — min-distance fields; iterate a 1-cell
   seam-strip exchange to convergence (order-independent, so distributed == serial). Same shape as the
   conduit-resolution boundary pass (§2), with `min(distance)` instead of a terminal label.
3. **Per-flat `max` reduction** for `flat_height[label]` — a cross-tile reduction keyed on the unified
   label.
4. **Boundary-aware mask + `d8_masked_FlowDir`** — a 1-cell mask halo suffices for the final local
   steepest-descent-on-mask (reuse richdem's `d8_flow_flats` / `d8_masked_FlowDir`).

### Why deferred (not a step too far — right home)
It is a distributed reimplementation of Barnes-2014 and reuses the MPI build's label-unification +
reductions; building it standalone in the in-process column-split harness would duplicate that
machinery and then discard it. Option 3 already bounds the problem, so this is safely additive
"in the middle, whenever."

### Acceptance criteria
- On any tiling (incl. tiles smaller than a giant lake), the flat branch is **bit-identical** to serial
  `flowdirs` (extends the current `FLOWDIR-MATCH` to the currently-capped cases), tree unchanged.
- Per-rank memory `O(N/P)+O(boundary)` (no growing halo); rounds ≈ flat diameter.
- Validated against the differential oracle exactly as the conduit pass was.

### Related
- Depends on: the resolve_flats flat-routing switch for serial DH (flag for Richard, plan §10).
- Bounded-for-now by: option 3 (`halo_cap`, DONE).
- Benchmarks + option-3 flag exposure (choosing cap / non-identical-but-functional modes): a follow-on
  task once the MPI build exists.

---

## ENH-2: seedless tiles have no exterior seed for PhaseAB

**Status:** open, characterized (found 2026-07-27 while validating the MPI harness increment 1).
Pre-existing — the in-process stitch (`tools/dephier_stitch.cpp`) has it too, not introduced by the
distributed build. Not a blocker for the seeded cases the core was validated on (fractals + edge
fixtures), but a robustness gap the "correct for ANY input" (option-3) philosophy wants closed.

**Type:** correctness / robustness (tiling seeding).

### Context
`GetDepressionHierarchyPhaseAB` throws `std::runtime_error("No OCEAN or BOUNDARY cells found, could
not make a DepressionHierarchy!")` when a tile contains **no ocean cell and no BOUNDARY cell**. The
stitch (and the MPI harness, identically) pre-labels a seam-edge cell BOUNDARY only if its steepest
descent *crosses the seam and is not to ocean*. A tile that is an interior closed bowl — every edge
cell drains inward, no ocean — therefore gets no exterior seed and PhaseAB aborts. Reproduced on
`kerry_test.dem` split 20→col 5 (2 ranks) and `kerry_test11.dem` split 7,15; both the stitch and
`dephier_mpi.exe` abort with the identical exception, while sibling fixtures of the same size/split
(e.g. `kerry_test1`, `kerry_test12`) have a seed and pass. So it is **data-dependent**, not geometric.
At 30″ global scale with large tiles a fully-seedless continental-interior tile is very unlikely (a
long perimeter almost surely has one outward-draining cell), but "very unlikely" is not "impossible",
and a build that must be correct for any input has to handle it.

### What it does / options
The physical fact: a closed-bowl tile still spills *somewhere* once filled; its lowest rim point is its
true outlet, and if that point is on a seam it belongs in the cross-tile outlet graph. Candidate fixes,
cheapest first:
1. **Open every non-true-grid-edge seam column as exterior when a tile would otherwise be seedless.**
   Seed the tile's seam-edge cells as BOUNDARY (provisional exterior) so PhaseAB has a seed; the conduit
   pass then resolves each to the real depression/ocean it drains into, exactly as for the normal
   cross-draining BOUNDARY cells. Mirrors parallel_priority_flood, where tile edges are *always* open
   boundaries, not conditionally. Likely the right, minimal fix.
2. **Guard + fallback:** detect a seedless tile and seed only its single lowest seam-edge cell as the
   provisional exit. Less uniform than (1); keep as fallback thinking only.

Either changes tile seeding, so it touches the **stitch too** and must be re-validated against serial
with the differential oracle (the tree must stay bit-identical). Flag alongside the other serial-output
-adjacent determinism items for Richard (plan §10).

### Acceptance criteria
- A seedless-bowl tile builds without throwing; its basin resolves to the same depression/ocean serial
  assigns (differential oracle STITCH-MATCH + FLOWDIR-MATCH), tree bit-identical.
- Regression fixture: a small DEM + split that currently aborts (e.g. `kerry_test.dem` split 5) becomes
  a passing `mpi_phaseab_*` / `stitch_edge_*` case.

### Related
- Shares the conduit-resolution machinery (§2) for resolving the newly-opened BOUNDARY cells.
- Found by: `tools/dephier_mpi.exe` increment 1 (per-rank PhaseAB), which faithfully reproduces the
  stitch abort — evidence the harness matches the stitch even in its failure modes.

---

## RESEARCH DIRECTION: lake ↔ drainage-network integration (CHONK-informed)

**Status:** parked (2026-07-27). A research direction, not a bounded code task. Build when ready;
intend to explore **joining forces with Boris Gailleton** (see below).

### The idea (Wickert)
The depression hierarchy (DH) and the tributary river network are both binary directed graphs, and
they *intersect* — rivers flow through lakes and wetlands. Splice them into one graph so depressions
act, on the network, as **sinks** (sediment), **capacitors** (water storage), and **mixers** (dissolved
matter incl. pollutants). Stop pretending the world is a single downhill-integrated network (or
digitally carving/filling DEMs to force it).

### Verification (what already exists — do NOT rebuild)
Checked the literature (searches + reading the actual paper/code). The splice + sediment routing are
**already done, single-node, and are our own lineage** — do not reinvent them:
- **CHONK 1.0** (Gailleton, Malatesta, Cordonnier, Braun; *GMD* 17, 71–90, 2024) builds its lake tree
  "with a principle adapted from Barnes et al. (2020)," cites Barnes 2021 (FSM) and Callaghan & Wickert
  (2019) as closest prior aim, and descends from Garcia-Castellanos' **TISC** (endorheic-LEM pioneer).
- The splice = a **"fake link"** from each base depression's pit cell to the cells just downstream of its
  outlet (Fig. 4a, §3.3.2), inserted so lakes sort before their downstream, then cancelled after the
  topological order is computed. `depressiontree.hpp`: `externode`/`internode`/`node2tree`/`linkhood` +
  level-sorted `get_treestack`.
- Sediment routing through lakes is solved via **"de-processing"** (§3.4.3): submerged cells' fluvial/
  hillslope fluxes are reversed and recomputed with the lake present, guarded by mass-balance checks —
  the authors call it "convoluted." A generic well-mixed-reservoir operator (`mix_two_proportions`,
  volume-weighted, renormalized) mixes **sediment provenance** proportions.
- Also covered: evaporative (leaky) capacity; multiple-flow-direction depressions. "Lakes in LEMs" is an
  active subfield (also Salles 2019; Camports 2020; Yuan 2019; Geurts 2018), not a gap.

### What is genuinely still open (our candidate contributions, on top of CHONK)
1. **Distributed / massively-parallel at global scale.** CHONK is single-node (ran on one 32 GB box) and
   "slow, lots of memory." Our distributed DH (tree+volumes+flowdirs, real MPI, footprint-bounded) is the
   differentiator — bring CHONK-style lake↔network integration to trillion-cell scale.
2. **Reactive dissolved solutes / pollutants.** CHONK has the mixing *primitive* (for particulate
   provenance) but tracks no dissolved chemistry — add concentration + a reaction/decay term per node
   (generalizes the per-lake nitrogen-removal limnology literature to the whole hierarchy).
3. **Richer lacustrine deposition.** CHONK's in-lake deposition is a simple uniform "draping," explicitly
   flagged by the authors as future work.
4. **Multiple flow directions (MFD)** — CHONK generalizes the DH to MFD; ours (like Barnes) is D8.

### Scope note
Our tool is *static-DEM analysis at scale*; CHONK is a *landscape-evolution model* (evolving topography).
Different animals — complementary, which is partly why a collaboration is natural.

### Collaboration
Consider joining forces with **Boris Gailleton** (Géosciences Rennes; also GFZ Potsdam). His CHONK →
**DAGGER** (graph backend) → **Scabbard** (topographic-analysis/LEM toolbox) stack is the natural partner,
and he is re-implementing CHONK's lake solver in that cleaner codebase. Our distributed DH + his
lake↔network machinery is a clean division of the four open items above.
- Flowdir interface note: raw D8 flow accumulation is *already solved and parallel* — richdem's
  `parallel_d8_accum` (Barnes 2017) consumes our flowdirs directly; don't rebuild it.
