# ENH-4: area-weighted volumes — port WTM's `cell_area` threading into the distributed DH

Draft enhancement, written like a GitHub issue (same register as `ENHANCEMENTS.md`). Held standalone to
stay clear of that file's in-flight parallelization revisions; fold in as ENH-4 (confirm the next free
number) at a good time. **Parked** by Wickert (2026-07-28): design is ready and the reference
implementation already exists — but it waits, being independent of the parallelization work in flight.

---

**Status:** parked / design ready (2026-07-28). **Not new work** — the area-weighted formulation is
already implemented and validated in **WTM** (our own lineage; see Provenance). This is a **port** of
WTM's `cell_area` interface into *this* repo's (distributed) DH, plus one widening (per-row → per-cell)
for the icosahedral case. Independent of ENH-1/2/3 and of the distributed stitch — it touches only the
Phase-D volume accumulation.

**Type:** feature / correctness (volumes on non-equal-area grids).

## Motivation
This repo's DH assumes **uniform unit cells**. Capacity is built from two *unweighted* accumulators:

```
include/dephier/dephier.hpp:884   cell_counts[clabel]++;                    // count of cells
include/dephier/dephier.hpp:885   total_elevations[clabel] += dem(i);       // unweighted Σ elevation
include/dephier/dephier.hpp:924   dep_vol = cell_count*out_elev - total_elevation;   // Σ_cells (out_elev−elev), area ≡ 1
```

Exactly right on an equal-area Cartesian grid, wrong on any grid where cell area varies — lat/lon (E–W
width ∝ cos φ), icosahedral/cubed-sphere (per-cell metric), or any projected/non-square raster. Upstream
richdem's FSM shares the same equal-area assumption (`DepressionVolume(sill, cell_count, total_elev)`;
`water_vol += wtd(c)`), so the gap runs through the whole surface-water stack.

## Provenance — WTM already solves this (do not reinvent)
`/home/awickert/models/WTM` forked DH + FSM and **added a per-cell area vector**, threading it through
exactly the accumulators above. It is the proven reference for this change:

- **Area is a first-class DH argument** — `GetDepressionHierarchy(topo, cell_area, label, …)`
  (`WTM/src/dephier.hpp:193,272`; called `WTM/src/WTM.cpp:237,471`, `run_dephier.cpp:100`).
- **Per-row area from latitude** — `cell_area[j] = cellsize_n_s_metres * cellsize_e_w_metres[j]`
  (`WTM/src/run_dephier.cpp:79`; `dvec cell_area` in `ArrayPack.hpp:63`).
- **Area-weighted volume** — `CalculateMarginalVolumes` (`WTM/src/dephier.hpp:752,765`):
  ```
  total_areas[clabel]   += cell_area[y];                          // area measure (= dep_area)
  total_volumes[clabel] += (out_elev - dem(x,y)) * cell_area[y];  // Σ (out_elev−elev)·area
  ```
  (tests carry the matching `area_times_elevation_total` accumulator.)
- **FSM water area-weighted to match** — `water_vol += runoff(x,y) * cell_area[y]`
  (`WTM/src/fill_spill_merge.hpp:327,413`), and groundwater storage/recharge/ocean-loss use `cell_area[y]`
  throughout (`transient_groundwater.cpp`).

So the equal-area limitation was **already fixed downstream in our own code**; upstream DH/FSM (and this
distributed fork) simply never inherited it. Same posture as the CHONK finding: port, don't rebuild.

## Machinery (the port)
The two DH forks have **diverged structurally** — WTM computes volumes in a separate
`CalculateMarginalVolumes`; this repo accumulates inline (884/885/924) and is the *distributed* build. So
this is a port of WTM's **pattern and interface**, not a wholesale lift:

1. **Thread an area accessor** through this repo's `GetDepressionHierarchy` (mirroring WTM's `cell_area`
   argument). Prefer an **accessor**, not a mandatory stored field, to protect the trillion-cell memory
   budget:
   - **lat/lon:** area is a function of *row only* → a per-row `cell_area[y]` vector, O(height) not O(N)
     (exactly WTM's `dvec`).
   - **icosahedral / general:** area varies per *cell* but is a deterministic function of grid position →
     **compute on the fly** from the cell metric, avoiding a full O(N) area array.
   - **default constant** → reduces to today's formula (see backward-compat).
2. **Area-weight this repo's accumulators** (the port of 804–805 onto 884/885/924):
   - add `dep_area += area(x,y)` (new; keep `cell_count` as the raw count),
   - `elev_area += elev·area(x,y)` (area-weighted analog of `total_elevation`),
   - `dep_vol = out_elev·dep_area − elev_area`.
   Touch-points: struct fields (`dephier.hpp:83–93`), leaf accumulation (884–885), tree aggregation
   (890–891, 916–919 — new fields sum over children identically), formula (924), and the `dep_vol`
   assertion (926).
3. **Distributed:** the area accessor is purely tile-local (each tile knows its rows'/cells' geometry);
   the area-weighted accumulators reduce over the tree exactly as the current ones do. **No new
   communication.**
4. **FSM (when/if it enters this repo):** WTM already shows the matching fix (`runoff·cell_area`); apply
   the same accessor there. Out of scope for the DH-only port but noted so it isn't re-derived.

## Acceptance criteria
- **Constant area ⇒ bit-identical** to current `dep_vol` (× the constant): differential oracle on the
  fractal/edge fixtures stays clean. (WTM's own tests pin this with `cell_area = {1,1,…}`.)
- **WTM as differential oracle for lat/lon:** on a shared DEM + per-row `cell_area`, this repo's ported
  volumes match WTM's `total_volumes` — a direct check against the proven implementation.
- **Per-cell (2D) area:** a synthetic where area varies *within* a row (icosahedral-like) matches an
  analytic volume — exercises the widening beyond WTM's per-row vector.
- Tree-aggregation invariant preserved (parent `dep_area`/`dep_vol` = Σ children).

## Related
- **Reference implementation: WTM** (`/home/awickert/models/WTM`, `src/dephier.hpp`,
  `src/fill_spill_merge.hpp`) — our lineage; this ports its `cell_area` interface and generalizes per-row
  → per-cell.
- **Area leg** of the geometry interface for the spherical / icosahedral work (spherical-grid tangent in
  `ENHANCEMENTS.md`); prerequisite for correct volumes on any non-equal-area global grid.
- **Equal-area grids sidestep it:** an exactly-equal-area icosahedral construction (e.g. Snyder) lets the
  accessor return a constant, reducing to today's formula — complementary to this port, and a point in
  favor of choosing an equal-area icosahedral grid.
- **Independent** of ENH-1 (flat resolution), ENH-2 (bowl-interior tiles), ENH-3 (internally drained):
  those touch seeding/flat/flow; this touches only Phase-D volumes. Schedulable entirely on its own.
