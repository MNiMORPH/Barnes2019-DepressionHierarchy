# Enhancements (deferred, designed) — distributed DepressionHierarchy

Local tracker for designed-but-deferred work on the parallel/distributed DH build, so the
plan stays on a linear path to completion. Each entry is written like a GitHub issue and can
be posted to the MNiMORPH fork's tracker when wanted.

---

## ENH-1: bit-identical O(boundary) distributed flat resolution ("option 2")

**Status:** DONE, now GENUINELY PER-RANK (2026-07-27 core `23bba02`/`0abce35`; 2026-07-29 distribution
`9075f99`). The label-free three-relaxation form (`resolve_flat_flowdirs_option2`, `tools/dh_flats.hpp`)
is bit-identical to richdem's `resolve_flat_flowdirs`, superseding option 3's halo cap, O(boundary) for
ANY flat extent. **The distribution plumbing is now built:** `resolve_flat_flowdirs_rank` runs the three
relaxations per rank on the tile + a 1-column halo, with a per-round seam exchange and a gather→OR→
broadcast convergence all-reduce (`FlatComm` in `dephier_mpi.cpp`); the central full-grid pass is retired.
No rank materializes a full-grid field. Validated: `flat_diff=0` (per-rank flat cells bit-identical to
serial) across all 335 fixture×split cases; suite 14/14; `tools/flat_mask_reconstruct_test.cpp` guards the
underlying identity. This was the last per-rank algorithm step; the whole distributed build (tree +
flowdirs) is now O(N/P)+O(boundary).

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

**Status:** **DONE (2026-07-28)** under the agreed acceptance bar — **correct volume + valid hierarchy**,
NOT byte-identical-to-serial (Andy's scoping call 2026-07-28). Every bowl-interior tile that used to abort
now builds pit-only with the right volume and a valid tree. Commits: `523b3f6` core flag, `7c358ce`
harness, `766a7bf` collapse Pass B, `263ab19` diagnostic guard, `261bbcd` NoData-as-ocean outlet fix,
`1fba9f9` `DH_AUDIT_OUTLETS`. At global 30″ this is a **guaranteed** case (large endorheic basins bigger
than a tile — Caspian, Tarim, Chad, Altiplano, Great Basin). The first candidate fix (open the seam as
exterior) was implemented and **reverted** (it severs tile-spanning basins); the correct base fix is
**pit-only seeding**; see Resolution.

### Implementation status vs serial (2026-07-28) — MEASURED, bar MET
Exhaustive fixture sweep `kerry_test*/testdem* × all split columns` at `-9999` (a crude uniform ocean level;
several fixtures then show *pre-existing* seeded DIFFERs unrelated to ENH-2 — a sweep artifact, ignore). Of
the **110** splits that used to ABORT (now build pit-only): **all 110 are volume-correct vs serial (0
mismatches) and produce a valid tree** (canonicalize throws on a cycle/multi-parent; it never does). So the
bar is met for every case. Byte-identity breakdown, for the record: 70 bit-identical (STITCH-MATCH); ~40
differ ONLY in `ocean_linked` nesting order — the documented **PhaseCD tie-break / out_cell** class (changes
serial output; already on Richard's review list), which the bar deems acceptable.

**MPI harness confirmed (2026-07-29):** the same fix went into `dephier_mpi.cpp` (`a8ecd2e`; its outlet
re-derivation is a separate copy — see ENH-5), and a **full MPI sweep** (`dephier_mpi.exe`, all 335
fixture×split cases) is **0 bad-volume, 0 error**: 73 MPI-TREE-MATCH, 262 DIFFER-but-volume-correct. So the
real distributed build meets the bar too — not just the in-process stitch.

**What landed and is validated (0 regressions vs the serial oracle on the full sweep):**
- `permit_without_baselevel_seed` flag + harnesses passing it (bowl-interior tiles build instead of abort).
- collapse `Pass B` meta-dissolve widening (+14 flat-straddle cases now bit-identical, incl. seeded ones).
- divergence-dump bounds guard (latent hang/segfault on stale post-collapse labels).
- **NoData-as-ocean outlet fix** (`261bbcd`): the re-derived outlet scan skipped all NoData cells, but
  `ocean_labels` makes NoData OCEAN, so basin→NoData-ocean outlets were dropped → `inf` volume (34 cases).
  **A PRE-EXISTING stitch bug** (reproduces on `ddb3eca`; minimal repro `testdem8.dem 0 1` → serial 239 /
  stitch inf), surfaced but not caused by ENH-2. Fixed; all 34 now volume-correct. Trade-off: 2 kerry_test4
  splits went MATCH→DIFFER, but audit-proven outlet-faithful and volume-correct — same tie-break class.
- **`DH_AUDIT_OUTLETS`** debug flag (`1fba9f9`): diffs PhaseAB's `tile.outlets` against the re-derived
  global set in real time, to catch future drift between the two outlet-discovery paths.

**Remaining (out of scope under the agreed bar; documented for completeness):** exact serial `ocean_linked`
nesting for the ~40 tie-affected cases = the PhaseCD tie-break/out_cell class. NOT a piecemeal fix — it
changes serial output and is best handled as one effort with Richard. Tracked, not blocking.

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

## ENH-3: permit DH in internally drained regions (no external base level)

**Status:** design sketched (2026-07-28, from the Mars / seafloor design exchange); **scoping DECIDED
2026-07-28 = Option A** (mode (b), closed-catalogue), shaped so a configurable base level (Option B) is a
purely *additive* extension later. Implementation pending and **gated on ENH-2** — this is the small
top-of-tree increment on ENH-2's pit-only seeding. Not a blocker for Earth-with-ocean.

**Type:** feature / generality (base-level model). A capability increase, not a robustness fix.

### Motivation
DH currently *requires* an external base level: PhaseAB seeds from OCEAN/BOUNDARY cells and PhaseC closes
the tree only through an OCEAN-endpoint outlet (`dephier.hpp:699`) — hard-wiring "the world drains to the
sea." The important internally drained problems have no sea at all:
- **Whole planets without an ocean.** Mars (MOLA) is one closed endorheic system; its ultimate low is
  Hellas (~−8 km). A global Mars run has zero ocean cells — it aborts today.
- **Earth's seafloor as the object of study.** Invert the usual base level: catalogue ocean-floor
  depressions (trenches, abyssal basins) rooted at the global bathymetric minimum (Challenger Deep,
  ~−11 km), instead of treating the sea surface as base level.
- **Isolated endorheic basins in their own frame.** A regional DEM clipped to the Caspian, Tarim, Chad,
  Altiplano, or Great Basin watershed is a *whole* closed domain, not a tile of a larger one.

The ocean-seeded case is the *special* case; internally drained is the general one.

**ENH-2 vs ENH-3 (precise):** ENH-2 = a **bowl-interior *tile*** inside an ocean-bearing domain (its rim
is in a neighbour tile, and the open-root chain lands at a real ocean elsewhere). ENH-3 = a **closed
*domain*** — no ocean anywhere, so that chain never terminates and the hierarchy must close on a virtual
base level. ENH-3 reuses all of ENH-2's pit-only flood; it adds only the domain-level closure.

### Resolution — Option A (mode b, DECIDED), with mode (a) as the quick alternative
**(b) Closed-catalogue — DECIDED.** Keep every basin, the deepest included, as a real depression. Reuse
ENH-2 `permit_without_baselevel_seed` so PhaseAB runs pit-only; the domain's top depression comes out
**open** (a dangling non-ocean root, per ENH-2's definition). The one genuinely new piece is the
**domain-level closure**: with no ocean anywhere, ENH-2's "chain neighbour → neighbour until it lands at
the ocean" never terminates, so the surviving open root is **adopted by the virtual ocean** (depression 0,
`pit_elev=−∞`, already created unconditionally at `dephier.hpp:354`) as a child at spill elevation **+∞**
— it never overflows. No real cell is dissolved; the ultimate low is simply the deepest *leaf*.

Distributed form: one global all-reduce — "does the domain contain **any** ocean cell?" If none, adopt the
dangling root(s) under ocean 0 at +∞ in the final PhaseC assembly. If the domain *does* contain ocean
(Earth), behaviour is byte-for-byte unchanged — mode (b) is inert.

**(a) Ocean-analog — zero code, documented alternative.** Relabel the single deepest cell (or a small
region) as OCEAN and run unchanged: the throw is satisfied and the flood roots there. Correct when the
ultimate low *is* the base level ("if Mars had a sea it would pool in Hellas"). **Caveat:** that cell
becomes base level, so the deepest basin dissolves into ocean-draining terrain — Hellas / Challenger Deep
is *erased* as a catalogued depression. Fine for an ocean-analog framing; wrong when the deepest basin is
the object (seafloor). Kept as the quick path; (b) is the general resolution.

### Why A, not the full base-level abstraction (Option B), now
Option B = make base level a first-class configurable spec `{ocean cells | ultimate-low | none/closed |
arbitrary datum z₀ | inverted}`. **Deferred, not rejected:** most of its extra reach is already free
without an API (mode (a) is a relabel; inverted/seafloor is running on −DEM), so over A it mostly buys
ergonomics, not capability — and it has **no concrete chosen-elevation consumer yet**, while costing more
churn on the published `dephier.hpp`. Decisive point: **A→B is additive, not rework**, *provided* A's
closure is written as "close against a base-level notion that today has two values (ocean cells |
virtual-ultimate)," not as a hardcoded no-ocean special-case. Build it that way and B slots in later (add a
datum case) with no back-tracking. Adopt B the moment a real chosen-elevation use appears (paleo-shoreline,
reservoir/fill level, a WTM datum coupling). (Global spherical grids — the earlier "Option C" — are out of
scope here; they belong to `GLOBAL_GRID_DESIGN.md`.)

**Code gate (trace before calling A "small"):** the +∞ adoption is a **post-outlet-loop parent
assignment** — there is no discovered OCEAN-endpoint outlet to trigger the `dephier.hpp:699` branch, so
closure is a *new terminal step* after the outlet loop, not a change to that branch. Confirm the insertion
trips no existing assertion. Two correctness items to verify, not assume:
- **N disconnected closed components**, not one: each connected closed component yields exactly one
  dangling root; adopt each under ocean 0.
- **Volume at +∞ must resolve to the finite "fill the component to its own max" volume**, not infinity —
  trace `CalculateMarginalVolumes` / `CalculateTotalVolumes`.

### Acceptance criteria
- A fully closed domain (no OCEAN/BOUNDARY — a single-bowl fractal, or a Mars-MOLA / seafloor subset)
  builds without throwing under mode (b); every basin including the deepest is catalogued; each closed
  component's tree has exactly one root (ocean 0) whose sole child is that component's top metadepression.
- **Serial** closed-domain DH (`permit_without_baselevel_seed` + top-closure) is the differential oracle:
  distributed == serial, tree + volumes bit-identical, exactly as for the ocean-seeded cases.
- Earth-with-ocean is byte-for-byte unchanged (the all-reduce makes mode (b) inert).
- Mode (a) documented as the quick ocean-analog path, with its deepest-basin-dissolves caveat stated.

### Related
- **Builds on ENH-2** (pit-only seeding, `permit_without_baselevel_seed`, the open-root definition).
  Implement ENH-2 first; ENH-3 is the top-of-tree increment on it.
- **Option B (configurable base level)** deferred pending a concrete chosen-elevation consumer; A is shaped
  so B is additive.
- Enables **seafloor depression studies** (root at the bathymetric minimum) — a distinct science use.
- The **planetary / global-sphere** case additionally needs `GLOBAL_GRID_DESIGN.md` (no boundary; for
  Mars, no ocean) — separate, later.

---

## ENH-5: consolidate the triplicated outlet re-derivation into one shared scan

**Status:** DONE 2026-07-29 (`b110f0d`). Extracted `tools/dh_outlets.hpp` (`OutletDB<CellIdx>` +
`outlet_skip` + `outlet_scan_intra`/`outlet_scan_seam`); the stitch, mpi oracle, and mpi distributed path
all call it. Pure refactor, validated bit-identical (stitch sweep MATCH 75 unchanged; mpi OUTLET/TREE/VOL
match; suite 22/22; real-MPI compiles). The consolidation itself SURFACED a drift bug it exists to prevent:
the stitch `record()` lacked the mpi's explicit lower-out_cell tie-break, making its outlet scan
order-dependent (intra-before-seam) — fixed first in `2948000` (MATCH 72→75; kerry_test4 splits 3/8/9
DIFFER→MATCH, previously **misattributed** to the Richard-coordinated tie-break class). The cross-rank
gather-merge `mrg()` keeps its own copy of the tie-break (a distinct reduction over already-reduced ORecs,
not a cell scan) — left as is. **Originally identified** 2026-07-28 while fixing the NoData-as-ocean inf
bug (ENH-2).

**Type:** refactor / correctness-hardening.

### Problem
The tiled build re-derives its outlet set from the resolved label grid (rather than reusing PhaseAB's
`tile.outlets` — see ENH-2's rationale for *why* it re-derives). But that scan is **triplicated**, three
independent copies that must agree with nothing enforcing it:
- `tools/dephier_stitch.cpp` — the `record()` intra + seam scans.
- `tools/dephier_mpi.cpp` **oracle** — `reduce(odb, …)` intra + seam.
- `tools/dephier_mpi.cpp` **distributed** — `reduce(myintra/myedgedb, …)` intra + seam.

The NoData-as-ocean inf bug was the same `if(isNoData) continue` in **all three**; the fix had to be
applied three times (`261bbcd` stitch, `a8ecd2e` mpi oracle + dist). Every future outlet-semantics change
is a 3× edit and a standing drift risk. `DH_AUDIT_OUTLETS` (`1fba9f9`) catches drift between PhaseAB and
the *stitch's* re-derivation at runtime, but it does not prevent the triplication.

### Fix
Extract the re-derivation into ONE shared helper (the way `dh_collapse.hpp` is shared), parameterized by a
label accessor, an elevation accessor, a `skip(cell)` predicate (NoData-but-not-OCEAN), and the
tile/seam iteration — so the stitch, the mpi oracle, and the mpi distributed path all call the same code.
Then a fix or a semantic change happens once. Validate against the differential oracle exactly as now
(0 regressions; full sweep + suite).

### Related
- Fixes the class the NoData bug (ENH-2) belonged to; all three already use the same `skip()` rule.
- `DH_AUDIT_OUTLETS` then guards the seam-exchange path rather than three parallel implementations.

---

## ENH-6: declared NoData semantics + base-level sill for OCEAN cells (eat real DEMs)

**Status:** open, characterized (2026-07-29, hit while building the Corsica/GEBCO example). Not a blocker
for the current examples (WTM marks its ocean with a sensible flat `0`, which is unambiguous); required to
consume real-world DEMs directly.

**Type:** correctness / robustness / generality (input handling). Touches shared core (`ocean_labels`,
outlet elevation) → coordinate with Richard (same review bucket as the tie-break items).

### Two distinct problems
**(a) OCEAN-cell elevation is taken from the raw DEM value (arithmetic fragility).** An outlet's elevation
is `max(landCell, oceanCell)`, and the code uses the ocean cell's *raw DEM value* for `oceanCell`. So the
DH's answer depends on the arbitrary NoData sentinel: sentinel `9` → outlet `9` (testdem8, "works" by
luck); `-9999` → outlet = the land value; **`NaN` → `max(land,NaN)=NaN` → PhaseC's outlet-sort assertion
fires and the build aborts** (observed on `corsica_gebco.tif`, GEBCO's native NoData is NaN). An OCEAN cell
is *base level*, so a land→ocean outlet's sill is the **land cell's** elevation; the ocean side should be
treated as −∞ (base level), never the sentinel. This makes the result robust to any sentinel, NaN included,
and independent of the sentinel's magnitude.

**(b) NoData is SILENTLY interpreted as OCEAN (semantic ambiguity).** `ocean_labels` does
`isNoData → OCEAN`. But on a real DEM NaN/NoData means several different things:
- **ocean / base level** (bathymetry masked out) — a terminal sink; ✅ the current assumption, sometimes.
- **interior void** (cloud, no-return lidar, a hole in the terrain) — NOT base level; treating it as ocean
  carves a false drainage path straight through the terrain.
- **outside the study area.**

The failure mode to kill: silently guessing "NoData = ocean" and then either crashing (NaN) or, worse,
**succeeding with a WRONG tree** (a weird sentinel quietly moving outlets).

### Fix
- **(a)** In the outlet computation (serial PhaseAB + the re-derivation scans + collapse), treat an OCEAN
  cell as base level: a land↔OCEAN sill = the land cell's elevation; never read the ocean cell's DEM value.
- **(b)** A declared-semantics flag, e.g. `--nodata {ocean | void | error}`:
  - `ocean` → base level (with (a)); `void` → an interior hole handled deliberately (fill / wall / carve —
    its own sub-question); default `error` → refuse to run if the DEM has NoData whose meaning wasn't
    declared. No silent interpretation.

### Acceptance criteria
- `dephier`/stitch/mpi build directly on `corsica_gebco.tif` (NaN NoData) under `--nodata ocean` with no
  sea-flatten pre-step, tree + volumes matching the flat-ocean-at-0 result.
- A DEM with an interior NoData void + `--nodata error` (or no flag) exits with a clear message, not a crash
  or a silently-wrong tree.
- Serial and all existing fixtures unchanged when NoData is absent or is a sensible ocean value.

### Update 2026-07-29 — a THIRD silent NoData-as-ocean read, in the collapse pass
Building the Corsica example split-invariant (identical tree regardless of tiling) uncovered another place
that silently reads NoData as OCEAN — the same class as problem (b), in a new spot. The §3.2 collapse pass's
`is_seam_artifact` (`tools/dh_collapse.hpp`) scans a degenerate pit's D8 neighbours for a cross-tile escape
but **skipped NoData neighbours** (`if(isNoData) continue`). richdem reads Corsica's ocean cells (value 0)
as `isNoData=1`, so a coastal degenerate pit whose only lower neighbour is the NoData-ocean across the seam
was never recognised as a seam artifact → left as a split-dependent extra zero-volume leaf (`pit(78,63)` at
split 78: 164 nodes vs serial 163). This is the identical `if(isNoData) continue` bug already fixed in the
two `record()` outlet re-derivation scans (`261bbcd`/`a8ecd2e`) — a **third copy** of the same NoData-ocean
assumption (see ENH-5: consolidate them; the whole family shares one `skip()`/"OCEAN = base level" rule).
Fixed `4c3edcb`: a cross-tile `isNoData || val<=pit_elev` neighbour is the escape; Corsica is now
bit-identical at every tiling. **Bearing on this ENH:** every one of these scattered reads bakes in "NoData
= ocean"; the `--nodata` flag (fix (b)) and the base-level sill (fix (a)) must land in ONE shared helper
(ENH-5) so a future `--nodata void` flips all of them at once, not three-plus places independently. The
collapse fix keeps the *current* `ocean_labels` semantics (NoData = ocean) — consistent, not a new policy.

Also added this session: the **node count is now an explicit decomposition-correctness diagnostic**
(`0371dec`, `STITCH-`/`MPI-DECOMP-CORRECT`/`-INCORRECT`). Because the tree is split-invariant, an unequal
node count is a NECESSARY-condition symptom that the split changed the tree. It flags TWO things (do NOT
read it as cleanly isolating decomposition bugs — an earlier framing that was wrong): (1) genuine seam
artifacts the collapse missed (spurious extra depressions — real bugs); (2) the STRUCTURAL sub-class of the
PhaseCD tie-break — at a TIED outlet, two basins rebuilt as a META vs one ocean_linked into the other (same
depressions + volume, but meta-vs-ocean_linked moves the node count; kerry_test2 4 vs 3). Only the pure-
ORDERING tie-break sub-class preserves node count (DIFFER but DECOMP-CORRECT, e.g. kerry_test4 3/8).
Progress this session on the split-invariance program the diagnostic measures:
- **Corsica coastal artifact fixed** (`4c3edcb`, the NoData-as-ocean escape above).
- **testdem8/kerry9 rim-fragment artifacts fixed** (`2e86ab7`, collapse Pass B2: dissolve a meta over {real
  basin leaf, one degenerate seam-artifact rim fragment} — the complement of Pass B's floor straddle).
- **Remaining (DECOMP-INCORRECT, ~35 fixtures, all volume-correct + valid-tree):** the meta-vs-ocean_linked
  tie-break class (class (2) above; appears at ALL splits, Richard-coordinated). A careful pass on the
  seam-dependent split-10 outliers (kerry10/11/12) found they are NOT a further collapse gap: the surviving
  seam artifacts are CELL-ASSIGNMENT differences (an artifact degen leaf carries cc=4 where serial's legit
  degen carries cc=1, and the parent meta's own aggregate is already wrong, e.g. cc 21 vs 18), so the
  tree-only collapse cannot repair them — it can dissolve nodes but not move cells between sibling basins.
  Same root as the tie-break class (the seam partitions a flat's cells differently than serial's flood at a
  tie); fix belongs upstream at the seam flat/tie partitioning, with the Richard-coordinated effort — NOT
  more collapse passes. The collapse (Passes A/B/B2) is complete for the tree-structural artifact classes.
  None is a volume regression.

### Related
- Surfaced by the Corsica example (README). The (a) sentinel-as-elevation is in shared core → Richard.
- The collapse NoData read (Update above) is the third instance of the (b) silent NoData-as-ocean class;
  fold all of them into ENH-5's single shared scan so `--nodata` governs them uniformly.

---

## ENH-7: geometry-deterministic flat LABEL assignment (a deterministic label set)

**Status:** RESOLVED IN-PROCESS 2026-07-29 — but NOT the way this entry first imagined (see below). The
"geometry-deterministic nearest-pit rule that changes serial labels" turned out to be the wrong model: the
flat label partition is ORDER-dependent (the flood's highest-index-first pop order), and the fix is to
REPRODUCE serial's exact flood partition, not replace it. Done via the flat-partition REPLAY
(`DH_FLAT_PARTITION_REPLAY` in `dephier_stitch.cpp`, proven partition-identical to serial on 49 DEMs by
`tools/flat_partition_replay_proof.cpp`), plus the one genuinely-serial-changing piece — the outlet
`out_cell` tie-break, now canonical (RICHARD_REVIEW_NOTES.md #3). Together: bit-identical serial↔tiled on
all 107 fixture×split cases (tag `serial-parallel-bitidentical`). The prototype's "nearest-pit, ~0 cost"
reading below was self-referential and is superseded. **Distributing this into the real per-rank MPI build
is ENH-8** (v1 rank-0 gather being built now; v2 = fully distributed). Original future-note kept below for
history.

**(historical) Status:** noted for the future 2026-07-29 (not planned); filed as **GitHub issue #3** to make it
permanent. Only needed if we ever want the distributed build's LABEL set (which depression each cell
belongs to) to be split-invariant — the current bar is volume-correct + valid-tree, which does NOT require
it.

**PROTOTYPE MEASURED (2026-07-29, `flat_relabel_prototype.exe`, commit `c263e61`) — the cost is ~ZERO, not
"changes serial output" as first assumed.** A geometry-deterministic rule (each flat cell → the nearest
depression pit over same-elevation D8 adjacency, tiebreak lower pit index) reproduces serial's flood labels
with **0 changes on every covered flat cell across all `test_cases`** (ocean 0). Serial's priority flood
already labels flat floors by nearest-pit, so the split-invariant geometric rule *is* what serial computes
— adopting it would barely touch serial. So ENH-7 may be nearly free, which strengthens the case. Two
caveats before believing it fully: (a) the tiebreak "lower pit index" coincided with serial's frontier
pop-order on these fixtures but is not *proven* identical universally — a flat with a genuine
equidistant-between-two-pits tie could differ by a handful; (b) "sill" flat cells (flat but above their
basin pit, or draining to ocean at equal elevation) are not covered by this simple pit-BFS and need
separate handling — though they are not the observed drift. Net: promising, worth a fuller pass, and no
longer clearly a "serial output changes" item — re-check the tiebreak on a tie-heavy fixture next.

**Type:** determinism / split-invariance of the per-cell label field.

### Problem
The depression a FLAT cell is assigned to is a byproduct of the flood's pop order (`dephier.hpp` labels a
cell as the first depression whose front reaches it). Richard's `radix_heap` fix (RICHARD_REVIEW_NOTES.md
#2) made the pop order deterministic *within a tile* (equal-elevation cells pop in cell-index order), so a
within-tile flat labels identically serial vs tiled. But **across a seam** the two floods still reach a
flat's cells in different global orders, so the tiled build labels some flat cells into a different
(equal-elevation, equally-valid) depression than the whole-grid flood. Measured residual (all VOL-MATCH,
valid trees): `kerry_test2` split 3 = 12 cells, `kerry_test10` = 209, `kerry_test11` = 407, `kerry_test12`
= 320. Where those relabelled cells sit on a saddle, they also re-binarize the tree (meta vs `ocean_linked`
at a tied outlet) — the DECOMP-INCORRECT-but-volume-correct residual (the "class 2" of the 2026-07-29
tie-break characterization; see RICHARD_REVIEW_NOTES.md and [[distributed-dh-stitch-state]]).

The SAME flood-order-on-flats root also drives the "class 3" cases (e.g. `kerry_test3`/`kerry_test7` split
3), where the labels match serial exactly (0 raw-label diffs) but the OUTLET CELL of a flat saddle is
chosen differently — serial's flood picks it by flood order, the geometric re-derivation by lowest index —
so `PhaseC` re-binarizes meta vs `ocean_linked` (verified 2026-07-29: serial 10 outlets vs re-derived 5,
on different elevation-6 rim cells). A geometry-deterministic flat rule that fixes the label assignment
would fix the outlet-cell selection too; both are the same problem.

### The idea (if we pursue it)
Replace flood-order flat labelling with a GEOMETRY-deterministic rule — the label analog of what
`resolve_flats` (Barnes-2014) already does for flat FLOWDIRS, and what ENH-1 distributes as
order-independent relaxations. A flat cell's owning depression would be a deterministic function of the
geometry (e.g. nearest lower outlet by the same away/towards gradient fields ENH-1 builds), independent of
which flood front arrived first. Then flat labels are split-invariant by construction, no cross-seam order
to reproduce. **Cost:** it CHANGES serial labels on flat-heavy DEMs (an improvement — removes latent
flood-order dependence — but Richard's to bless, same bucket as the radix fix). Distributes over the same
seam-exchange machinery ENH-1 uses. Only worth it if a downstream consumer needs a stable per-cell label
field (e.g. exact lake-basin masks that must not shift with the tiling); the tree + volumes are already
correct without it.

**NOTE (superseded, 2026-07-29):** this "geometry-deterministic, changes-serial" idea is NOT what we did.
Serial's flat label is order-DEPENDENT (highest-index-first flood pop order), so we REPRODUCE it exactly
(the replay) rather than replace it. See the status block above and ENH-8.

---

## ENH-8: distribute the flat-label REPLAY per-rank ("v2", fully O(cap·boundary))

**Status:** designed 2026-07-29. **v1 (rank-0 gather-resolve-scatter) is being built now** in `dephier_mpi`;
this entry is **v2 = the fully-distributed, no-rank-0-bottleneck form**, filed to keep the plan on a linear
path. Not a blocker: v1 already gives the correct (serial-identical) distributed result; v2 is the footprint
optimization, exactly the way `PhaseCD` is v1-centralized with the fully-distributed 2016 join deferred.

**Type:** performance / footprint (v1 is correctness-complete).

### Context
The flat-label partition fix (ENH-7, resolved in-process) is the **replay**: reproduce serial's exact flood
label partition over a connected flat, seeds + the radix (elevation, index) pop order. Proven distributable
as a **per-tile adaptive-halo replay** (`tools/flat_label_distributed.cpp`, DIST-PARTITION-MATCH on the whole
corpus + size-200 fractals up to 5 tiles). Two keys: (1) label each basin by its PIT's GLOBAL cell index —
namespace-free, so straddling basins agree across tiles with **no seam stitch / union-find**; (2) convergence
needs TWO consecutive stable halo doublings (one is a coincidence on a wide flat).

**Why this needs a WIDE (growing) halo, unlike ENH-1's O(boundary) flat FLOWDIRS.** ENH-1 distributes as
three ORDER-INDEPENDENT relaxations (1-column seam exchange per round, converge regardless of sweep order).
The flat LABEL partition is ORDER-DEPENDENT — it *is* the flood's sequential pop order (a flat interior cell
becomes its own pit if it pops before a wavefront reaches it). So it cannot be independent relaxation rounds;
it is an ordered replay over the flat's extent, needing a halo ≈ the flat's cross-seam width (O(flat extent),
the same bound as the flowdir option-3 halo). For fixtures with flats wider than a tile, that halo spans
MULTIPLE tiles.

### v1 (being built): rank-0 gather-resolve-scatter
Mirror the existing conduit pass (`dephier_mpi.cpp` ~line 438, "v1 gather-and-resolve"): each rank ships its
seam-band columns (dem + glab) to rank 0; rank 0 runs the windowed pit-index replay over each seam-straddling
band and scatters back corrected labels; the existing outlet scan / PhaseCD run on the corrected labels.
Footprint: O(band·boundary) gathered to rank 0 (not the whole grid). Correctness-complete, codebase-consistent.

### v2 (this enhancement): fully per-rank wide halo
Each rank fetches its OWN adaptive/capped halo (up to the flat extent) via a **multi-column seam exchange
chained across tiles** (a halo wider than one tile pulls columns from the neighbour's neighbour, etc.), runs
its own pit-index replay, and relabels its owned flat cells to `glab` at the pit cell (exchanged alongside
the dem halo). No rank-0 gather. Footprint per rank O(N/P)+O(cap·boundary); cap the halo (option-3 style) and
accept a valid-but-maybe-not-identical partition for flats wider than the cap — exactly the flowdir tail.

### Acceptance criteria
- `MPI-TREE-MATCH` on the corpus at multiple tilings, per-rank memory O(N/P)+O(cap·boundary), no whole-grid
  or wide-band gather to rank 0.
- Bit-identical to the v1 result within the cap; documented fallback beyond it.

### Related
- Supersedes/continues **ENH-7** (the flat-label determinism, resolved in-process via the replay).
- Sibling of **ENH-1** (the flat-FLOWDIR distribution) — but ENH-1 is O(boundary) relaxations; ENH-8 is
  O(cap·boundary) ordered replay, because labels are order-dependent and flowdirs are not.
- Proven algorithm: `tools/flat_label_distributed.cpp`; feasibility: `tools/flat_seam_feasibility.cpp`.
- Serial side already done: the `out_cell` tie-break (RICHARD_REVIEW_NOTES.md #3) + the replay
  (`serial-parallel-bitidentical` tag).

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
