# WTM API review — the DH public interface, from an integrator's outside eyes

External review of this repo's public depression-hierarchy API, written from the perspective of the
**WTM** consumer entering from outside. The framing question: *could the local DH become a drop-in
replacement for the richdem-lineage DH that WTM currently vendors, and where is the API unfamiliar or
improvable to a caller who already knows the richdem original?*

Written by the WTM-side integrator (Claude, with A. Wickert), 2026-07-29. **Review/opinion only — no
code changed here.** Companion to `../RICHARD_REVIEW_NOTES.md`, `PARALLEL_DEPHIER_DESIGN.md`, and the
already-parked `../ENH-4-area-weighted-volumes.md`.

## TL;DR

The **PhaseAB / PhaseC / PhaseCD decomposition is the good idea** and is well-drawn — it puts the seam
in the right place (the tree is the one inherently-global object; only Phase C must be serial). But the
serial wrapper is **not a drop-in** for the richdem/WTM signature today: five API/type divergences stand
in the way, four of which are choices that would make the API *more familiar to any richdem consumer*,
not just to WTM. The highest-leverage move is to give the convenience wrapper the **exact richdem
signature as a superset**, keep the phased API as the new power-user surface, and let WTM eventually
**delete its vendored `../src/dephier.hpp` and consume this header** — converging the two forks instead of
maintaining both.

## Context: WTM vendors its own DH

WTM does **not** consume the richdem submodule for depression hierarchy. It carries a modified copy at
`WTM/src/dephier.hpp` (namespace `richdem::dephier`) plus `WTM/src/fill_spill_merge.hpp`. FillSpillMerge
(FSM) is WTM-side and is **not** in this repo — so any integration swaps **only the DH build** and keeps
WTM's FSM untouched. The one interface that must line up is the `Depression` struct that
`GetDepressionHierarchy` produces and FSM consumes.

## What the new API gets right

- **The phase split.** `GetDepressionHierarchyPhaseAB` (grid-reading flood + outlet discovery,
  per-tile parallel) / `PhaseC` (grid-*free* tree assembly, inherently central) / Phase D
  (`CalculateMarginalVolumes` + `CalculateTotalVolumes`, distributable) is the correct decomposition.
  The tree is the unavoidably-global object and it is isolated to the one phase that has to be serial.
- **Deterministic outlet sort** (`PhaseC`, `out_elev → out_cell → pit-cell minmax`) is the load-bearing
  idea: it makes a distributed build reproduce the serial tree bit-for-bit at triple junctions. This is
  *exactly* the determinism the WTM serial≡parallel equilibrium demo depends on — the two projects are
  aligned on what matters most.
- **`permit_without_baselevel_seed`** for bowl-interior tiles is a thoughtful distributed affordance.
- Keeping a one-call `GetDepressionHierarchy` wrapper over the phases is the right courtesy to serial
  callers.

## Why it is not a drop-in today — five gaps

| # | WTM (richdem lineage) | Local DH | Kind |
|---|---|---|---|
| 1 | `GetDepressionHierarchy(dem, `**`cell_area`**`, label, `**`final_label`**`, flowdirs)` | `(dem, label, flowdirs)` | signature |
| 2 | `dh_label_t = int32_t` | `uint32_t` | type (see below) |
| 3 | area-weighted volumes via `cell_area`; `dep_area` field | uniform-cell `dep_vol = cell_count·out_elev − total_elevation` | physics — **WTM addition, = your ENH-4**, not a local regression |
| 4 | emits `final_label` (highest metadepression per cell); FSM consumes it | not produced | output surface |
| 5 | struct carries `wtd_vol, wtd_only, my_cells, dep_area` (FSM/GW storage) | carries `total_elevation` instead | struct |

**On #3 — the honest attribution.** Stock/richdem DH (and therefore this repo) is **cell-weighted** and
correct on an equal-area grid. The area-weighting is a **WTM-side enhancement** (`WTM/src/dephier.hpp`
threads a per-cell `cell_area` vector through the volume accumulators), already captured here as the
parked `../ENH-4-area-weighted-volumes.md`. So this is not a local regression; it is a known port. It still
blocks a *drop-in*, because consuming the local wrapper as-is on a WTM lat/lon grid (E–W width ∝ cos φ)
yields wrong capacities — the budget will not close. Landing ENH-4 closes it.

## The phase letters do not match the paper

Separate from the WTM drop-in, but surfaced while reviewing the API: **the code's phase names do not
match the algorithm's own paper.** Barnes, Callaghan & Wickert (2020, *Earth Surf. Dynam.* 8, 431–445),
§3 "The algorithm", defines four phases:

- **A — §3.1 Ocean identification** (seed the PQ with ocean cells)
- **B — §3.2 Pit cell identification** (seed the PQ with pit cells)
- **C — §3.3 Depression assignment** (priority-flood: label cells to leaf depressions; record the lowest
  outlet between each adjacent depression pair)
- **D — §3.4 Hierarchy construction** (sort outlets by elevation; union-find sweep low→high to merge
  leaves into meta-depressions)

Volumes are **not** one of the four phases — they are **§6.4 "Depression statistics"**, explicitly
post-processing.

The code merged the paper's A+B+C into one function and shifted the letters by a whole phase:

| Paper §3 | Code today | What the name implies vs. what it is |
|---|---|---|
| **A** Ocean identification | ⎫ | |
| **B** Pit-cell identification | ⎬ all inside `GetDepressionHierarchyPhaseAB` | `PhaseAB` reads as "ocean + pit seeding" (paper A+B) but also contains the entire flood + outlet discovery (paper C) |
| **C** Depression assignment | ⎭ | |
| **D** Hierarchy construction | `GetDepressionHierarchyPhaseC` | **called C, is the paper's D** |
| §6.4 statistics (volumes) | `...PhaseCD` extra step / `Calculate{Marginal,Total}Volumes` | **called D, is not a lettered phase at all** |

The **decomposition is right** — the parallelization seam sits exactly where it should: grid-touching
work is paper A–C (+ marginal volumes, distributable); grid-free work is paper D (+ total volumes,
central). That C|D boundary is precisely the code's `AB` | `C` split. Only the **labels** are wrong, and
they have already caused ambiguity in the design docs ("does `PhaseC` mean phase C or phase D?").

**Recommendation:** rename to the paper's terms and drop the letter-soup — e.g.
`FloodAndAssignDepressions` (paper A–C, one grid pass), `ConstructHierarchy` (paper D, grid-free),
`DepressionStatistics` / `CalculateVolumes` (§6.4). If letters are kept, make them the paper's A/B/C/D
and note that the implementation fuses A–C into a single pass. A paper-faithful rename is the natural
moment to also give the convenience wrapper the richdem signature (rec. 1) — one API-cleanup pass.

## Recommendations, ranked (to make it drop-in / more familiar)

1. **Give the serial wrapper the exact richdem signature, as a superset.**
   `GetDepressionHierarchy<elev_t,topo>(dem, cell_area, label, final_label, flowdirs)`. Keep the current
   phased API as the added power-user surface; let the convenience wrapper mirror the interface everyone
   already knows. Closes #1 and #4 at once and makes it a recognizable drop-in. (Concrete signature in
   the appendix.)
2. **Revert `dh_label_t` to `int32_t`.** `uint32_t` is the most surprising change to a richdem reader,
   and it is a latent bug: `LastLayer` still assigns `-3` (`../include/dephier/dephier.hpp:946`), which on a
   `uint32_t` wraps to ~4.29e9. Signed labels + a named sentinel is the familiar and correct form.
3. **Land ENH-4 (area-weighting).** For WTM this is physics, not familiarity — it feeds the
   conservative-flux, machine-zero water budget on non-equal-area grids. Aligning its signature to the
   original `cell_area` argument (rec. 1) *is* the delivery vehicle for ENH-4.
4. **Decide where the FSM/GW fields live** (`wtd_vol, wtd_only, my_cells, dep_area`). Cleanest long-term:
   the local `Depression` is a **superset** carrying them — zero cost to a pure-DH consumer that ignores
   them — so WTM can delete its vendored copy and consume this header. One `dephier.hpp`, not two drifting
   ones.

Minor / consistency nits: `Depression` changed `class → struct` (cosmetic); `total_elevation` is a new
public field worth a units/intent line; the `-3` magic in `LastLayer` wants a named constant regardless
of #2.

## The 2nd-order Picard+Anderson coupling

The WTM implicit groundwater solve never touches the DH directly. The coupling is: DH+FSM sets the
lake / `wtd_vol` configuration each cycle, and the 2nd-order-in-time implicit solve (semi-implicit
Picard / BDF2-on-V, with Anderson as the matrix-free default) equilibrates against it. What the implicit
solve *requires* of the DH is that its output be **rank-deterministic** — otherwise WTM's serial≡parallel
property collapses. The local DH's deterministic outlet sort delivers exactly that, so the design is
aligned. The one caveat is #3/ENH-4: the volumes must be area-weighted, or the discrete budget will not
close on a real lat/lon grid.

## Strategic recommendation

Do **not** reconcile by editing WTM's vendored copy. Push recs 1–4 into *this* repo so the serial
wrapper is a superset-compatible drop-in for the richdem signature. Then WTM's integration reduces to
"delete `../src/dephier.hpp`, point the include at this header, leave `fill_spill_merge.hpp` untouched," and
the two headers **converge** instead of forking further. That serves the parallel-DH goal and leaves one
source of truth.

## Appendix — the superset wrapper signature

```cpp
// Convenience wrapper: richdem-compatible superset. Wraps PhaseAB -> PhaseCD, threads cell_area
// into Phase D (ENH-4), and writes final_label in the marginal pass. The 3-arg phased API stays
// for power users / distributed builds.
template <class elev_t, Topology topo>
DepressionHierarchy<elev_t> GetDepressionHierarchy(
    const Array2D<elev_t>&     dem,
    const std::vector<double>& cell_area,     // per-row (or per-cell) area; unit vector reproduces stock behavior
    Array2D<dh_label_t>&       label,         // dh_label_t == int32_t
    Array2D<dh_label_t>&       final_label,   // output: highest metadepression per cell
    Array2D<int8_t>&           flowdirs);
```

Passing a unit `cell_area` reproduces today's cell-weighted result exactly, so the superset is
backward-compatible with the equal-area path.

## Verification

Every assertion above was checked against source on 2026-07-29 (not recalled). Line references:

**Local DH** (`../include/dephier/dephier.hpp`):
- `dh_label_t = uint32_t` — line 31 (gap #2).
- `LastLayer` assigns `mylabel = -3` where `mylabel` is deduced `uint32_t` (from `label(x,y)`), then
  writes it back to `label` — lines 940, 946, 950. On `uint32_t` this stores ~4.29e9 (rec. 2). A leftover
  from the signed-label era, now inconsistent with the uint32 switch.
- `struct Depression` (not `class`) — line 45. Fields present: `cell_count` (83), `dep_vol` (87),
  `water_vol` (90), `total_elevation` (93). **Absent**: `wtd_vol`, `wtd_only`, `my_cells`, `dep_area`
  (gap #5).
- Cell-weighted volume `dep_vol = cell_count·out_elev − total_elevation` — line 924; accumulators
  `cell_counts[clabel]++` / `total_elevations[clabel] += dem(i)` — lines 884–885 (gap #3, stock behavior).
- No `final_label` anywhere; `CalculateMarginalVolumes(deps, dem, label)` takes no `cell_area` —
  lines 856–860 (gaps #1, #4).
- Phase API: `PhaseAB` (284, `permit_without_baselevel_seed` at 290), `PhaseC` (645, deterministic sort
  `out_elev → out_cell → pit-cell minmax` at 660–667), `PhaseCD` (817), 3-arg serial wrapper (841–851).
- `BOUNDARY = max−1` distributed-tile sentinel — line 193 (not a gap; noted for context).

**WTM vendored DH/FSM** (`WTM/src/`):
- `dh_label_t = int32_t` — `dephier.hpp:29`.
- Area-weighting is WTM-added: `cell_area` threaded through `CalculateMarginalVolumes`
  (`dephier.hpp:763–822`, `total_areas`/`dep_area`/area-weighted `total_volumes`) and `CalculateTotalVolumes`
  (`842, 847, 851`); called with `(topo, cell_area, label, final_label, flowdirs)` at `WTM.cpp:247,501`.
- `final_label` is **consumed**, not just produced: written in the marginal pass (`dephier.hpp:794`) and
  read at `fill_spill_merge.hpp:536` inside `CalculateWtdVol` (`fill_spill_merge.hpp:202`), the routine that
  builds `wtd_vol`/`wtd_only`.
- ENH-4 attribution: `../ENH-4-area-weighted-volumes.md` ("port of WTM's `cell_area` interface into this
  repo", parked 2026-07-28) — the area-weighting is a WTM enhancement, not a local regression (gap #3).

**Paper** (phase-naming section): Barnes, R., Callaghan, K. L., and Wickert, A. D.: "Computing water flow
through complex landscapes – Part 2: Finding hierarchies in depressions and morphological segmentations",
*Earth Surf. Dynam.* 8, 431–445, 2020 (https://doi.org/10.5194/esurf-8-431-2020). §3 headings verified
verbatim against the published HTML: 3.1 Ocean identification / 3.2 Pit cell identification /
3.3 Depression assignment / 3.4 Hierarchy construction; volumes at 6.4 Depression statistics.
