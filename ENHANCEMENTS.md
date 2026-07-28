# Enhancements (deferred, designed) — distributed DepressionHierarchy

Local tracker for designed-but-deferred work on the parallel/distributed DH build, so the
plan stays on a linear path to completion. Each entry is written like a GitHub issue and can
be posted to the MNiMORPH fork's tracker when wanted.

---

## ENH-1: bit-identical O(boundary) distributed flat resolution ("option 2")

**Status:** DONE (2026-07-27, commits `23bba02` core, `0abce35` harness integration). The label-free
three-relaxation form (`resolve_flat_flowdirs_option2`, `tools/dh_flats.hpp`) is bit-identical to
richdem's `resolve_flat_flowdirs` on the edge fixtures AND all kerry fractals, and is now the harness
flat pass (MPI-FLOWDIR-MATCH), superseding option 3's halo cap. It is O(boundary) for ANY flat extent.
*Remaining (same class of gap as the per-cell MPI verification):* the harness runs the relaxations on
the whole grid in-process; a real-MPI build runs the identical D8-stencil relaxations per tile with a
1-column seam exchange per round (structurally evident, not yet a separate driver). See the VALIDATED
SIMPLIFICATION note below for the mechanism, and `tools/flat_mask_reconstruct_test.cpp` for the proof.

Original framing (kept for context): designed, deferred (the problem was already *bounded* by option 3;
this recovers bit-identity on the giant-flat tail without growing per-rank memory).

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

### VALIDATED SIMPLIFICATION (2026-07-27) — no union-find, three relaxations
Reading `resolve_flats_barnes` closely, `label_this` floods equal-elevation D8-neighbours, so **"same flat
label" ⟺ "D8-adjacent + same elevation."** The two gradient BFSs already only cross same-label cells, so the
label check can be replaced by an **elevation check** — and `flat_height` (max away-dist per flat) becomes a
**max-relaxation over same-elevation adjacency**, also label-free. So the whole thing is **three
order-independent relaxations, no global flat labeling** (item 1's union-find is unnecessary):
  1. `away_dist` — multi-source min-BFS from high edges, over {NO_FLOW, same-elevation} adjacency.
  2. `towards_dist` — multi-source min-BFS from low edges, same adjacency.
  3. `flat_height` — max-relaxation of `away_dist` over same-elevation adjacency.
Then `mask = flat_height − away + 2·towards` and `d8_masked_FlowDir(mask)`, both purely local.
Each relaxation distributes as: local compute → 1-column seam exchange → relax → iterate to global
convergence (all-reduce of a "changed" flag); rounds ≈ flat diameter in tiles, O(boundary)/round.
**Proven bit-identical to richdem's `flat_mask`** on a full-grid reconstruction test (scratch
`flat_reconstruct.cpp`): `mask_diff=0` and label-free `flat_height diff=0` (302 flat cells, 65 high / 5 low
edges). Correctness foundation done; remaining ENH-1 work is the distribution plumbing (per-rank
relaxation + seam exchange + convergence loop, edge classification with the 1-column halo, integration into
`resolve_flat_flowdirs`'s capped branch, validation against the full-grid form).

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

## ENH-2: bowl-interior tiles have no base-level seed for PhaseAB

**Status:** design converged (2026-07-28), implementation pending. Found 2026-07-27 while validating
the MPI harness increment 1. Pre-existing — the in-process stitch (`tools/dephier_stitch.cpp`) has it
too, not introduced by the distributed build. Not a blocker for the seeded cases the core was validated
on (fractals + edge fixtures), but a robustness gap the "correct for ANY input" philosophy wants closed —
and at global 30″ a **guaranteed** case, not a rare one: large endorheic basins much bigger than a tile
contain fully bowl-interior tiles (Caspian, Tarim, Chad, Altiplano, Great Basin). The first candidate fix below (open the seam
as exterior) was implemented and **reverted** — it severs tile-spanning basins. The correct fix is
**pit-only seeding**; see Resolution.

**Type:** correctness / robustness (tiling seeding).

### Context
`GetDepressionHierarchyPhaseAB` throws `std::runtime_error("No OCEAN or BOUNDARY cells found, could
not make a DepressionHierarchy!")` when a tile contains **no ocean cell and no BOUNDARY cell**. The
stitch (and the MPI harness, identically) pre-labels a seam-edge cell BOUNDARY only if its steepest
descent *crosses the seam and is not to ocean*. A tile that is a **bowl interior** — every edge
cell drains inward, no ocean — therefore gets no base-level seed and PhaseAB aborts. Reproduced on
`kerry_test.dem` split 20→col 5 (2 ranks) and `kerry_test11.dem` split 7,15; both the stitch and
`dephier_mpi.exe` abort with the identical exception, while sibling fixtures of the same size/split
(e.g. `kerry_test1`, `kerry_test12`) have a base-level seed and pass. So *which* tiles are bowl interiors
is **data-dependent** (basin locations × tiling), but that *some* occur at 30″ global scale is
**guaranteed**, not rare: a large compact endorheic basin much bigger than a tile contains tiles whose
*entire* perimeter drains inward (Caspian, Tarim, Chad, Altiplano, Great Basin). The "long perimeter
almost surely has one outward-draining cell" intuition holds only when the basin is comparable to or
smaller than the tile; a basin bigger than the tile defeats it — exactly the regime of the big closed
basins the WTM cares about.

### The key realization: a "seed" does two separate jobs
PhaseAB's priority queue is seeded by two kinds of cell that play *different* roles — the serial case
hands you both together, so they look like one thing:

- **Pit seeds** — interior local minima (a cell with no lower neighbour). They **start the flood**; each
  becomes a *leaf* depression. Always present: any closed region has a bottom.
- **Base-level seeds** — the exterior (ocean/BOUNDARY) cells on the flood front. They **anchor the tree to
  base level**: where the flood's exterior region meets a land depression, outlet discovery emits an outlet
  with an OCEAN endpoint, and that is the only thing letting PhaseC close the tree at the root
  (`dephier.hpp:699`, the `depa_set==OCEAN || depb_set==OCEAN` branch).

A **bowl-interior tile** lies entirely inside a closed bowl whose rim is beyond the tile on every side:
every seam cell drains *inward*, so nothing is labelled BOUNDARY and there is no ocean. It has the
**bottom but not the rim** — pit seeds (the flood can run) but no base-level seed (no local anchor). The old
`exterior_cells==0` throw conflated the two jobs: it aborts for lack of the *anchor* even though the
*flood* needs only the pits.

### Resolution: pit-only seeding + let the existing seam machinery supply the outlet
Relax the need for a *local* base-level outlet. When the caller signals a bowl-interior tile (via
`permit_without_baselevel_seed`), PhaseAB runs **pit-only**: skip the throw, flood from the interior pits, and let the tile's top
depression come out **open** — no outlet to the exterior found within the tile.

"Open" is a property of the **outlet graph**, not a flag on the depression object (every leaf leaves
PhaseAB with `out_cell`/`odep`/`parent` unset — PhaseC fills them). Precisely:

> A depression is **open** ⟺ following `parent` up to the root of its connected component, that root is
> **not** ocean-linked (no ancestor has `ocean_parent`) and is not node 0 — a **dangling, non-ocean
> root**. A closed hierarchy has exactly one root: OCEAN. Note `ocean_parent==false` on a single node is
> the *norm* (a nested pond spills into its parent meta, not the sea); openness is the absence of an ocean
> link **along the whole path to the root**, i.e. at the component root.

The tightened statement of how it closes:

> Locally we allow a tile's top depression to be **open** — no local outlet to the exterior, so under
> local-only assembly its parent stays empty and it dangles as a non-ocean root. A **neighbour tile then
> supplies the missing seam edge** (the rim saddle, found by HandleEdge from the 1-column halo we already
> exchange). The **global PhaseC assembly uses that edge to give it a parent**, chaining neighbour →
> neighbour until the chain lands at the ocean. The neighbours supply the *edge*; the central assembly
> assigns the *parent* (a meta-depression born in assembly, not necessarily a depression in the adjacent
> tile).

So the whole change is: **PhaseAB doesn't need the rim to do its job (find the bottom + internal
structure); only the assembly needs the rim, and the assembly already looks across seams.** Sill
*elevations* are absolute DEM values from outlet discovery — the base-level seed affects **topology** (what
anchors to the root), not **geometry**.

**Scope of the code change (small):** add a `bool permit_without_baselevel_seed=false` parameter to
`GetDepressionHierarchyPhaseAB`; the default `false` preserves the serial throw exactly. Replace the
`exterior_cells==0` throw (`dephier.hpp:350`) with pit-only seeding when `permit_without_baselevel_seed`. Depression 0
(ocean) is already added unconditionally (`dephier.hpp:354–363`). The distributed harness/stitch drop the
bowl-interior abort and pass `permit_without_baselevel_seed=true`. Everything downstream — namespace remap
(a bowl-interior tile's empty ocean-0 is simply not shipped; the single global ocean is anchored from rank
0), HandleEdge, PhaseC,
PhaseD, collapse — is untouched.

#### Rejected (implemented, reverted): open the seam as exterior
Seeding a bowl-interior tile's seam cells as BOUNDARY (`dh_seed.hpp`, now deleted) **severs tile-spanning
basins**: the seam floor becomes a terminating sink, so the bowl-interior tile drains "out" at a spuriously low
elevation *and* neighbour cells draining across resolve into that sink. On `kerry_test.dem` split 5 the
tree shape came out right but `total_dep_vol` was serial=80 vs dist=5 — the volume above the false base
destroyed. Seen from PhaseC, it injected **false OCEAN-endpoint outlets** at the seam floor, so assembly
anchored the basin to a false low root. Pit-only injects *no* false exterior, so PhaseC anchors only where
a real base-level seed exists. Kept visible as the lesson: the failure was exactly at the seed→PhaseC anchoring
touchpoint, which is why "open the seam" and "pit-only" look similar but are opposites.

Either approach changes tile seeding, so it touches the **stitch too** and must be re-validated against
serial with the differential oracle (the tree must stay bit-identical). Flag alongside the other
serial-output-adjacent determinism items for Richard (plan §10).

### Acceptance criteria
- A bowl-interior tile builds without throwing; its basin resolves to the same depression/ocean serial
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

---

## PARKED TANGENT: DH on spherical / geographic grids

**Status:** parked (2026-07-27, Wickert). For later. Relevant to the *actual* global-30″ goal — the
current core is Cartesian, but a true global run is on a lat–lon (geographic) grid.

Three things a spherical build must handle:
1. **Periodic longitude.** The grid wraps at ±180°; the east and west edges are neighbours. Our seam/
   tiling machinery already resolves vertical seams — a global wrap is essentially "one more seam" joining
   column 0 to column W−1, so this may be a modest extension of the existing conduit/outlet/HandleEdge
   passes rather than new machinery.
2. **Latitude-varying cell geometry.** E–W cell width scales with cos(latitude); N–S ≈ constant. This
   changes (a) cell **area** → Phase D marginal **volumes** (`dx·dy` per row, not constant — feeds the
   already-distributed volume pass) and (b) any distance-weighted flow metric (D8 diagonal-vs-cardinal
   weighting). The tree/labels are unaffected; volumes and metric-weighted routing are.
3. **Pole singularity.** Meridians converge; the top/bottom rows are near-singular — the genuinely hard
   part (cap, exclude, or a pole-aware topology).

The tree topology is largely geometry-agnostic; the real work is volumes (per-row area), a periodic-x
seam, and poles. Worth scoping against real GEBCO 30″ (which is geographic) before a production run.
