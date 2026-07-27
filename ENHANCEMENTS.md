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
