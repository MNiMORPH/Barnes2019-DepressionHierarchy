# Distributed DepressionHierarchy — engineering plan for the MPI build

**Status (2026-07-27).** The in-process research core is complete and footprint-bounded, validated
bit-identical to serial (1521/1521), tagged `distributed-dephier-core`. This document lays out the
remaining work: turning the validated in-process algorithm into a real MPI build. It is the concrete
expansion of `PARALLEL_DEPHIER_PLAN.md` §9 step 5 ("MPI harness"), which that plan treats as one line.

Everything here is *engineering*, not research: every novel algorithm is already proven in
`tools/dephier_stitch.cpp` against the canonical oracle. The risk now is in the transport layer and the
distributed protocols, not in whether the stitch is correct.

---

## 1. Component decomposition

The single "MPI harness" line is really nine components. In dependency order, with status:

| # | Component | Status | Notes |
|---|-----------|--------|-------|
| 1 | **Transport layer** (`TiledArray2D` or reuse Barnes' `A2Array2D`) | **DECISION PENDING** (§3) | tile bounds, neighbour-rank map, per-cell channels (elev/label/flowdir), edge-strip extract+exchange. Everything sits on this. |
| 2 | **Per-rank Phase A/B** | reuse | run `GetDepressionHierarchyPhaseAB` on the local tile; edit = feed tile edges (BOUNDARY pre-label) not grid edges. |
| 3 | **Perimeter-strip exchange** | new, small | ship 1-cell `(elev,label)` strips to neighbour ranks (Barnes' HandleEdge model). |
| 4 | **Cross-tile outlets** | logic done, wire | HandleEdge over exchanged strips → each rank's contribution to the outlet DB. In-process code exists in the stitch. |
| 5 | **Conduit resolution, distributed** | **design below (§2)** | the boundary-graph pass, realized across ranks. New protocol; postdates the plan's §9. |
| 6 | **Global namespace remap** | designed (plan §4) | `Allgather` per-tile depression counts → prefix-sum offsets → remap every label field. |
| 7 | **Centralized Phase C** | designed (plan §2), **gated by §8** | gather outlet set + depression records to rank 0, run serial `PhaseCD` verbatim. |
| 8 | **Collapse + Phase D** | collapse done; D partly | collapse runs on the assembled tree (rank 0, O(#deps)); marginal volumes distributed, total reduced. |
| 9 | **MPI oracle harness** | pattern done | run distributed over localhost ranks, gather, diff canonical signature vs serial — the same bit-identity check used throughout. |

---

## 2. Distributed conduit resolution — MPI protocol

The in-process pass (`tools/dephier_stitch.cpp`, "conduit resolution") is two phases: (1) per-tile-local
`localwalk` reduces each BOUNDARY cell to a terminal or a seam exit; (2) `chase` follows exits across
tiles until a terminal. Measured shape (instrumentation, since reverted): BOUNDARY cells O(boundary)
(~300–2800 per seam), conduit paths ≤ 6 tile-crossings even at 8 tiles. This is a **shallow boundary
graph** on O(boundary) nodes — the property the distributed protocol exploits.

### Setup (already implied by components 2–3)
Each rank `r` owns tile `r`, has run PhaseAB (local labels, flowdirs, BOUNDARY seeds), and has exchanged
its 1-cell perimeter strips `(elev, local_label)` with neighbour ranks. Global label offsets from the
prefix sum (component 6) are known, so a local label maps to a global one.

### Phase 1 — local, no communication
Each rank runs `localwalk` for each of its BOUNDARY cells, using only its own `label`/`fd` arrays plus
the exchanged strip (for the `drain()` cross at a seam seed). Each BOUNDARY cell resolves to one of:

- **`RESOLVED(gid)`** — the local walk reached a local depression or ocean; `gid` is its global label.
- **`EXIT(cell_g)`** — the walk left the tile at a seam; `cell_g` is the *global cell id* of the entry
  cell just across the seam (in a neighbour's edge strip, so its coordinates and owning rank are known).

This is the exact `LW` struct the in-process Phase 1 produces, but keyed on a global cell id.

### Phase 2 — resolve the boundary graph
The entry cell of an `EXIT` is either a non-BOUNDARY cell in the neighbour (its label is a terminal,
readable from the exchanged edge strip) or a BOUNDARY cell whose own status is `RESOLVED`/`EXIT`. So the
graph is: BOUNDARY cells as nodes, `EXIT` edges pointing to the next cell, terminals as sinks. Two
realizations, chosen by the same trigger as centralized-vs-distributed Phase C (§8):

- **v1 — gather-and-resolve (recommended, pairs with centralized Phase C).** Each rank ships its
  O(boundary) `{global_cell_id → RESOLVED gid | EXIT cell_g}` records to rank 0 via `MPI_Gatherv`,
  alongside the outlet set it already gathers for Phase C. Rank 0 runs the in-process `chase` over the
  merged records (a pointer-follow to a terminal; shallow, so trivially fast), producing the final
  `global_label` for every BOUNDARY cell, and scatters them back. This reuses the centralized machinery,
  is deterministic by construction, and the gathered boundary graph is O(total boundary) ≪ N — the same
  footprint class as the outlet set already being centralized. **This is the natural v1.**

- **v2 — iterative propagation (only if §8 says the central gather is too big).** Fully distributed
  pointer-jumping: each rank publishes, for its edge BOUNDARY cells, their current resolution; each round
  a rank follows one hop for its unresolved `EXIT`s using neighbours' published edge resolutions; an
  `EXIT` whose target is now `RESOLVED` becomes `RESOLVED`. Chains of length ≤ k resolve after k rounds;
  since max chain (`maxcross`) ≤ 6, it converges in a handful of rounds. Terminate on an `MPI_Allreduce`
  of "any change this round" == 0. This is the analogue of Barnes' 2016 boundary-graph ("master graph")
  step, applied to conduits. Communication is O(boundary) per round, nearest-neighbour, ≤ ~6 rounds.

### Determinism (why either realization is bit-identical to serial)
The resolution is a pure function of the flowdirs and the `drain()` tie-break (lowest downhill neighbour,
ties by highest cell index) — both already deterministic and identical between serial and tiled builds
(that identity is what the 1521/1521 validation established). The gather-and-resolve realization runs the
*same* `chase` code the in-process harness runs, so it reproduces the in-process result exactly; the
iterative realization computes the same fixed point. Either way the output `global_label` grid is
byte-for-byte what the validated in-process stitch produces → the canonical signature matches serial.

### Message format
Per BOUNDARY cell: `(uint64 global_cell_id, uint8 kind, uint64 payload)` where `payload` is the global
label (`kind=RESOLVED`) or the entry global_cell_id (`kind=EXIT`). O(boundary) records per rank; a single
`Gatherv` (v1) or nearest-neighbour exchange (v2). **Global cell id (decision, Wickert): a single global
row-major index `y*W + x`, widened to `uint64`.** Rationale: (a) a single integer is simplest in practice
(indexing, sorting, hashing); (b) it is **bit-identical to the serial build's own cell indexing** —
`GetDepressionHierarchy` already uses the global flat index, so keeping it is what lets the distributed
result match serial cell-for-cell, as the in-process harness already does; (c) `uint64` removes the
`uint32` ceiling the plan (§4) flags at 30″ (only ~4.6× headroom; `uint64` has ~2×10¹⁰). A
`(tile_id, local_flat)` pair is a useful *debugging view* but not the on-wire identity — it would not
match serial's indices, forfeiting (b). (Resolves plan open-question §10.5.)

---

## 3. Transport-layer foundation — reuse vs. build (RESOLVED: reuse the *program pattern*, not the tile class)

A code map of Barnes' richdem tiling (`A2Array2D.hpp`, `Layoutfile.hpp`,
`programs/parallel_priority_flood/main.cpp`) found a clean separation that decides this:

- **`A2Array2D` is a single-process, out-of-core tile cache** — LRU eviction (`A2Array2D.hpp:91`),
  tiles paged to/from **disk** via a `Layoutfile` (`_LoadTile`, `dumpData()`, lines 104–151); no
  `mpi.h`, no ranks. It is templated but **single-channel** (`std::vector<std::vector<WrappedArray2D>>`,
  line 89), so carrying elev + label + flowdir needs **three parallel instances**. It is therefore
  **not** the MPI distribution abstraction — it is a within-process disk-cache helper.
- **richdem's MPI distribution lives in the `parallel_priority_flood` *program*** (Producer–Consumer:
  rank 0 orchestrates and holds the centralized `mastergraph`; ranks ≥ 1 each flood one tile;
  `CommInit`/`CommISend`/`CommRecv`). A2Array2D is not even used there — tiles are tracked by a manual
  `TileGrid`.
- **`HandleEdge`/`HandleCorner` already cross elevation *and* labels** (`main.cpp:344–398`), and
  **`label_offset`** already makes per-tile label namespaces global by prefix sum (`main.cpp:431–438`)
  — exactly our §4 remap.
- **`mastergraph`** (`main.cpp:421`) is a per-global-label spillover graph gathered on rank 0 and
  resolved in a single pass (lines 501–543) — structurally identical to our centralized Phase C +
  gather-and-resolve conduit design.

**Verdict / efficient integration path.** The efficient move is to reuse the *program pattern and the
edge/label machinery*, not the tile class:

1. **Do not build a heavyweight MPI-native `TiledArray2D`.** Follow the `parallel_priority_flood`
   Producer–Consumer choreography — it is Barnes' proven MPI structure, and our v1 (gather the outlet set
   + boundary graph to rank 0, resolve, scatter) is the direct analogue of his `mastergraph`. His code
   already validates the centralized-gather approach at scale for filling, which de-risks our centralized
   Phase C.
2. **Reuse `HandleEdge`/`HandleCorner` + `label_offset`** for components 3/4/6 — they already carry the
   depression-label channel we need across seams and globalize the namespace.
3. **`TiledArray2D` shrinks to a thin per-rank wrapper** (tile bounds, edge-strip extract, global↔local
   index) around richdem `Array2D` — not a distributed abstraction. Use `A2Array2D` per rank *only if* a
   single rank's tile exceeds RAM (out-of-core within a rank); a later refinement, not v1.
4. **What is genuinely new** (beyond his filling code) is only the richer centralized payload —
   depression *records* + the outlet *database* (vs. his mastergraph of spill elevations), the collapse
   pass, and the conduit boundary graph (§2). Everything else is his pattern with a heavier message.

Net: the integration is *efficient* — we inherit the MPI plumbing, edge exchange, and namespace
machinery, and only extend the centralized resolution payload. This is worth confirming with Richard
(plan §10.3), but the code evidence is unambiguous that his program pattern — not A2Array2D — is the base.

---

## 4. Real-tile footprint (§8) — does centralized Phase C fit? **VERDICT: yes, with 2–3 orders of margin.**

Measured on real GEBCO 30″ tiles (converted `gebco_08.nc` → GeoTIFF; contrasting terrains cut with
`gdal_translate -srcwin`; sea `z≤0` and the 1-cell border ring opened to ocean so each tile has a base
level; `tools/dephier_stats` on the unmodified serial `GetDepressionHierarchy`). Depression *density* =
depressions per land cell:

| Tile (15°-wide unless noted) | land cells | depressions | density | tree bytes/land-cell |
|---|---|---|---|---|
| Andes (high relief)          | 2,371,478 |  31,896 | 1.34% | 1.3 |
| Great Basin (endorheic, 9°)  |   772,838 |   9,933 | 1.29% | 1.2 |
| Australia interior           | 2,151,416 |  43,950 | 2.04% | 2.0 |
| Central US                   | 2,154,004 |  52,280 | 2.43% | 2.3 |
| Central Asia (30°×20°)       | 7,948,839 | 194,231 | 2.44% | 2.3 |
| Fennoscandia (glacial lakes) | 1,634,439 |  45,943 | 2.81% | 2.7 |
| Tibet/Himalaya               | 1,736,878 |  54,384 | 3.13% | 3.0 |
| Sahara                       | 2,154,004 |  81,393 | 3.78% | 3.6 |
| Canadian Shield (lakes)      | 1,731,136 |  78,050 | 4.51% | 4.3 |
| W. Siberia (flat/wet)        | 2,154,003 | 103,164 | 4.79% | 4.6 |
| Amazon (low relief)          | 2,153,997 | 128,478 | 5.96% | 5.7 |
| **aggregate (11 tiles)**     | **26,963,032** | **823,702** | **3.05%** | **~2.9** |

**Key findings:**
- Density is **1.3–6.0% of land cells**, land-weighted mean **~3%**; the outliers are low-relief/flat/wet
  terrain (Amazon, W. Siberia, Shield), not the high-relief or endorheic cases.
- Density is a **stable terrain property, not a function of tile size**: the 30°×20° Central-Asia tile
  (2.44%) sits squarely among the 15° tiles, so extrapolating per-cell density to the globe is valid.
- `Depression` node = **96 bytes**.

**Global projection** (N = 933,120,000 cells; land ≈ 29% ≈ 272M cells):
- depressions ≈ 272M × 3% ≈ **8.3M** (mean), or ≈ **16M** at a pessimistic all-Amazon 6%.
- light-tree footprint = deps × 96 B ≈ **0.8 GB** (mean) to **1.5 GB** (worst-case density).
- outlet database ≈ a few × deps × 16 B ≈ **0.4–1 GB** (plus transient build-time hash overhead).
- boundary graph (conduit records) = O(boundary) ≪ #deps — negligible.
- **Centralized Phase C working set ≈ 1–3 GB globally.**

**Verdict.** The object that must be centralized (the light tree + outlet set + boundary graph) is
**~1–3 GB** — 2–3 orders of magnitude below the O(N) full-grid footprint (DEM+label+flowdir+PQ ≈ 8–15 GB
*per the whole grid*) that forced distribution in the first place. So the split holds cleanly: **distribute
A+B (the memory wall, one tile per rank, O(N/P)); centralize C on rank 0 (fits on any modern node with
huge margin).** The v1 architecture — centralized Phase C + gather-and-resolve conduits (§2 v1) — is
confirmed viable; no need for the distributed-C / v2-conduit escalation.

*Honest caveats (do not change the order of magnitude):* the open-border tiles create some edge basins a
seamless global build wouldn't (slight over-count); marking `z≤0` as ocean drops the handful of
below-sea endorheic basins (Caspian etc.) (slight under-count). Both are small relative to the 2–3 orders
of headroom. A confirming measurement on the *actual* global build (once it exists) is worth doing, but
the architecture decision is not close.

---

## 5. Readiness verdict & recommended order

The three gates I flagged before starting are now largely settled:

1. **§8 footprint gate — PASSED (§4).** Centralized Phase C is ~1–3 GB globally, 2–3 orders below the
   distributed O(N) grid. v1 architecture confirmed; no distributed-C escalation needed.
2. **Transport-layer foundation — RESOLVED (§3).** Reuse Barnes' `parallel_priority_flood` Producer–
   Consumer program pattern + `HandleEdge`/`label_offset`; `TiledArray2D` shrinks to a thin per-rank
   wrapper. A2Array2D is a single-process disk cache, not the MPI layer. *(Worth a one-line confirm with
   Richard per plan §10.3, but the code evidence is unambiguous.)*
3. **Conduit protocol — WRITTEN (§2).** v1 (gather-and-resolve, mirroring his `mastergraph`) chosen; v2
   iterative propagation held in reserve, not needed given the §8 result.

**So: the design is ready to jump in.** The one genuine external dependency is a brief co-scope with
Richard (integration reuse, plan §10) — worth a message, but not a blocker to starting the harness.

**Recommended build order** (each step diffed against serial via the oracle, as the in-process core was):

1. **Thin `TiledArray2D` per-rank wrapper** + `Layoutfile`-style tiling metadata (component 1).
2. **2-rank localhost MPI**: per-rank PhaseAB → HandleEdge strip exchange → gather outlets + depression
   records + boundary-graph records to rank 0 → remap → serial PhaseCD → collapse → scatter/verify
   (components 2–9, v1). Diff canonical signature vs serial (plan §7 step 3).
3. **Scale ranks/tiles**; then **distributed Phase D volumes** (plan §9 step 6).
4. Escalate any single component to its distributed variant only if a later, larger footprint
   measurement on the real global build demands it (it won't, per §4).

---

## 6. Per-cell flowdirs across seams — status + the flat determinism case

`GetDepressionHierarchy` returns a per-cell `flowdirs` array alongside the tree. The distributed build
must reproduce it too (FSM routes runoff along flowdirs into the depressions).

**Non-flat crossings — DONE (bit-identical).** A tile flood cannot point a cell across its own boundary,
so seam-crossing cells kept a tile-local direction. The stitch's flowdir fix-up (commit `c47cb9e`)
restores serial's choice from the conduit pass's cross-seam data: a land pit/seam-seed points to its
lowest cross-tile neighbour ≤ the pit (highest-index tie); a sea-draining cell points to its highest-index
adjacent ocean cell. Result: **flowdirs are bit-identical to serial across the entire fractal sweep**
(β 1.3–2.5, seeds 1–15, single/multi-seam) — 0 diffs wherever flow direction is physically determined.

**Flat plateaus crossing a seam — KNOWN GAP (the determinism case).** On a large flat, DH assigns each
cell's direction as an **order-dependent byproduct of the flood** (`dephier.hpp:524` — the cell points to
whichever neighbour claimed it, set by the global radix pop order). When a seam cuts a flat, each tile
runs that claim from only its share of the flat's exits, so the in-flat routing differs — and the
difference propagates across the whole flat, not just the seam (measured 177/6084 cells on an 80×80
flat-plateau fixture; **0 on continuous terrain**). The **tree and every cell's sink are bit-identical**;
only the arbitrary in-flat direction (a gradient-less surface) differs. So this is a determinism/
completeness gap, not a hydrology error — FSM fills the same depressions with the same water either way.
It surfaces only on large flats (lakes, quantized plains) — the §8 vast-flats regime.

**Tractable path (Wickert): adopt Barnes-2014 `resolve_flats` instead of the flood byproduct.** richdem
already ships it (`richdem/flats/flat_resolution.hpp`, `Barnes2014.hpp`): `BuildAwayGradient` (BFS from
the flat's high edges) + `BuildTowardsCombinedGradient` (BFS toward its low edges) build a `flat_mask`,
and `d8_masked_FlowDir` reads each cell's direction off that mask. This routing is a **deterministic
function of the flat's geometry**, independent of processing order — so it is *reproducible*, which is
exactly what the flood byproduct is not. Cost to make flat flowdirs bit-identical distributed:
  1. **Switch DH's flat flowdirs to `resolve_flats`** (serial + tiled). Changes serial's flat directions
     (from arbitrary flood-claim to Barnes-2014 convergent routing — an *improvement*, and the standard);
     the depression **tree is unaffected** (flat directions don't touch labels/outlets). This is a core
     change, so flag for Richard alongside the other serial-output-changing determinism fixes.
  2. **Distribute the two gradient BFSs** with a cross-seam boundary exchange — the *same* pattern as
     conduit resolution (§2): each tile computes local high/low-edge distances, exchanges seam values,
     iterates to convergence (rounds ~ flat diameter in tiles; O(boundary) per round). Then
     `d8_masked_FlowDir` is a purely local read of the mask.
This is far cheaper than the intractable alternative (reproducing the flood's global pop order across
seams). Caveats: mesas (flats with no low edge) — `resolve_flats` handles them; a plateau spanning many
tiles needs more iteration rounds, but each is footprint-bounded. **Recommendation:** adopt
`resolve_flats` for flat routing (it improves the serial output too) and distribute it via the boundary-
exchange BFS when flat-flowdir bit-identity is wanted; until then, the gap is documented and harmless
(tree + sinks exact).
