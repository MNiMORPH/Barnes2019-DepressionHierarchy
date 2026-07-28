# Global (spherical) grid design for distributed DepressionHierarchy

**Status:** parked research direction (2026-07-28). Design/scoping notes, not an implementation plan —
per the Cartesian-first sequencing (§8), most detail is deliberately deferred; this document's job is to
preserve the *reasoning and decisions* and frame the *open questions* so they survive across sessions.

**Supersedes** the lat/lon "PARKED TANGENT: DH on spherical / geographic grids" section in
`ENHANCEMENTS.md` (its pole/periodic-longitude/area content is folded into §2, §5, §7 here). That section
should be **removed** from `ENHANCEMENTS.md` once the in-flight parallelization revisions settle, so the
spherical material lives only here.

Provenance of decisions: grid-basis preference is **Richard's** (icosahedral); the merits analysis and
sequencing were worked out with **Wickert**. From-memory external references (climate-model grid names)
are tagged for confirmation.

---

## 1. Purpose & status
The real global-30″ goal needs a global grid; the current core is Cartesian. This is the parked design
space for taking distributed DH to a whole planet (Earth *and* Mars). It is not started yet — it opens
when the Cartesian build is complete (§8), because finishing Cartesian is what maps the geometry interface
the port depends on.

## 2. Why not lat/lon — "spherical ≠ lat/lon"
The obvious "spherical grid" is a lat/lon (geographic) mesh, and it carries two problems that are
**artifacts of that particular tessellation**, not of doing DH on a sphere:
- **Pole singularity.** Meridians converge; the top/bottom rows are near-singular — a whole collapsing
  row, the genuinely hard case.
- **Latitude-varying area.** E–W cell width ∝ cos(latitude); this feeds Phase-D volumes (see ENH-4) and
  any distance-weighted routing. Tree topology is unaffected — geometry enters only at volumes/metric.
- **Periodic longitude.** The grid wraps at ±180° (col 0 adjacent to col W−1) — essentially "one more
  seam," a modest extension of the existing seam machinery rather than new work.

Choosing a *different* tessellation removes the first two by construction. That is the point of §3.

## 3. Grid-basis comparison & decision
Criteria that matter for a **depression/flow** method specifically:
- **(a) cell uniformity** — volume fidelity, less shape-induced bias in what counts as a depression;
- **(b) singularity mildness** — how bad, how many, how contained;
- **(c) flow-routing isotropy** — our entire output is directions + nesting, so azimuthal bias is a
  first-order quality concern.

| Grid | Area uniformity | Singularities | Stencil / isotropy | Notes |
|------|-----------------|---------------|--------------------|-------|
| lat/lon | poor (cos φ) | 2 singular *rows* (poles) | quad D8 (cardinal/diag bias) | native data grid |
| cubed sphere | ~2:1 (gnomonic; better equiangular) | 12 edge *kinks* + 8 three-way corners | quad D8 (inherits flat anisotropy) | **reuses Cartesian core + D8 ecosystem** |
| **icosahedral / hex** | near-uniform (~1.1–1.4 optimized) | 12 isolated pentagon *points*, no singular lines | hex, 6 equidistant nbrs (isotropic) | different connectivity; no D8 ecosystem |
| HEALPix | *exactly* equal-area | isolatitude | twisted D8 at face boundaries | equal-area ⇒ volumes need no weighting |
| reduced Gaussian | mitigated lat/lon | poles | isolatitude | not an escape from lat/lon |
| Yin–Yang | overlap needs interpolation | — | — | **ruled out** (interpolation breaks mass/topology exactness) |

**Decision: icosahedral on the merits (Richard's preference), cubed sphere as the fallback.** Setting
aside implementation ease (which favors the cube — see below), the hexagonal icosahedral grid wins on all
three criteria: most uniform area; only isolated pentagon point-defects (no kinked seams or collapsing
rows); and 6 equidistant neighbours at uniform 60° spacing remove the D8 cardinal-vs-diagonal length bias
that a flow-direction method is most sensitive to. The 60°-vs-45° angular-resolution tradeoff is real
(hex quantizes to 6 directions), but the literature read is that uniform neighbour distance beats
finer-but-anisotropic spacing for flow bias — **confirm with a source before committing to writing.**

**Where the cubed sphere still legitimately wins** (the fallback case): it stays quad-**D8**, so it reuses
our Cartesian per-tile core *and* the entire D8 hydrology ecosystem — notably richdem's `parallel_d8_accum`
(Barnes 2017), which consumes our flowdirs directly. Going hex makes that connectivity different
everywhere (flowdir encoding, accumulation). That interoperability is a genuine cost of icosahedral, not
just an ease argument — but it does not outweigh (a)+(b)+(c) for a depression/flow method.

## 4. Decomposition / MPI — does icosahedral tile cleanly? (Yes.)
The reassuring result: on the decomposition axis icosahedral is ~on par with the cubed sphere, via the
**rhombic layout**:
- Pair the 20 icosahedral triangles into **10 rhombic "diamond" panels**; each diamond is a **logically
  rectangular** i,j block after bisection. A rhombic (60°) lattice gives each interior cell **6 nearest
  neighbours at fixed index offsets** (±e₁, ±e₂, ±(e₁−e₂)) — clean hex connectivity inside a *rectangular
  memory array*. Each diamond subdivides into sub-blocks for more ranks.
- This is a solved problem at extreme scale — NICAM / ICON / MPAS decompose icosahedral grids this way
  (*from memory; confirm the specific model attributions*).

So the tile primitive (rectangular block + 1-cell halo) and the seam-exchange survive. **Two deltas** from
flat-Cartesian tiling:
1. **Tile adjacency is a fixed irregular graph**, not a 2D lattice; some diamond seams join with a
   transposed/reversed index orientation → a small set of precomputed **per-seam index transforms**.
2. **12 pentagon points** (5 diamonds meet; 5 neighbours) → bounded O(1) corner special cases.

Both are **also true of the cubed sphere** (12 edge flips + 8 three-way corners), so decomposability does
**not** favor the cube. If anything hex load-balances better: uniform cell area ⇒ uniform work/tile,
whereas the cube's ~2:1 area spread does not.

Our specific passes port: union-find flat labelling, conduit/outlet resolution, and the O(boundary)
relaxation-to-convergence all operate on *abstract* seam-adjacency, and their order-independence (what
makes distributed == serial) is untouched — they need only a well-defined seam adjacency, which the index
transforms supply.

## 5. The three geometry-interface legs
Geometry enters the DH in exactly three places; the tree topology is geometry-agnostic.
- **(a) Area / volumes — ENH-4, largely solved.** Area-weighted volumes are a **port of WTM's `cell_area`
  threading** (`/home/awickert/models/WTM`), generalizing per-row (lat/lon) → per-cell (icosahedral). See
  `ENH-4-area-weighted-volumes.md`. Prefer an area *accessor* (per-row vector for lat/lon; analytic from
  the metric for icosahedral) over an O(N) stored field, to protect the trillion-cell memory budget.
- **(b) Neighbour stencil / flowdir encoding — deferred.** Hex needs a 6-direction flowdir encoding
  (rhombic offsets); richdem's `parallel_d8_accum` will not consume hex flowdirs, so accumulation needs a
  hex analog or a mapping.
- **(c) Singularities — deferred.** Poles (lat/lon), pentagons (icosahedral), or corners (cube).

**Discipline (do-it-once without over-abstracting):** keep the Cartesian core **concrete** (D8, 2D-lattice
seams) through the correctness-critical ENH-1/2/3 work — a mid-stream pluggable stencil risks
destabilizing exactly the bit-identity we are proving. Instead produce a **"geometry interface map"** as a
byproduct of finishing Cartesian: the enumerated list of every area/stencil/singularity touch-point, so
the port is a bounded refactor, not a rewrite. Bake in now only the cheap, clearly-right leg — the
per-cell area accessor (leg a).

## 6. Scale & data
**Sizes (whole planet, best available global products):**

| | Product | Physical res | Cells N | Memory |
|---|---|---|---|---|
| Earth | GEBCO 15″ | ~463 m | ~3.7×10⁹ | ~75 GB |
| Mars | MOLA / MOLA–HRSC blend, native | ~200 m | ~5.6×10⁹ | ~110 GB |

Both are **single-fat-node feasible**; distribution buys speed, commodity-node portability, and the
future 1″ path (~8.4×10¹¹ cells ≈ 17 TB uniform; far less with sparse ocean tiles for land studies).

**Complexity / scaling:** DH is ~O(N log N) (O(N) with the linear-time priority-flood variant); it
distributes with **near-ideal weak scaling** (fixed cells/rank → fixed wall-time), with a communication
overhead growing only ~√N (surface-to-volume). At TB scale the **real wall-clock bottleneck is I/O and
DEM-mosaic assembly**, not the hierarchy computation.

**Data caveat:** GEBCO and MOLA are both **lat/lon** — a non-lat/lon target grid requires **regridding**,
which for a depression method is delicate (interpolation can create/erase pits/flats). An **equal-area**
target grid also simplifies volumes (leg a reduces to a constant). See §7.

**Internally drained (ENH-3):** Mars has *no ocean* — one closed endorheic system rooted at Hellas
(~−8 km); Earth's seafloor as an object roots at Challenger Deep (~−11 km); a global closed sphere has no
boundary at all. Required for any Mars run. See `ENHANCEMENTS.md` ENH-3.

## 7. Open questions (deferred until the interface is real)
1. **Data regridding** — conservative, depression-preserving resampling from lat/lon sources onto the
   target grid.
2. **Hex flowdir encoding + accumulation** — 6-direction encoding; the `parallel_d8_accum` ecosystem gap.
3. **Cell metric** — per-cell area (leg a) and neighbour distances; pentagon cells as exceptions.
4. **Global addressing** — (panel, i, j) vs a space-filling-curve linear index for comm/cache locality.
5. **Seam/halo bookkeeping** — per-seam orientation transforms; how a 1-cell halo behaves at a 5-way
   pentagon junction.
6. **Grid generation** — construction choice (recursive bisection / spring-optimized / Snyder
   equal-area); subdivision level matching ~15″ / ~200 m; the **equal-area vs shape-uniformity tradeoff**
   (an exactly-equal-area grid sidesteps most of leg a).
7. **Output & interop** — storing/visualizing on a non-raster grid, or reprojecting back to lat/lon.
8. **Validation** — analytic sphere test surfaces with known depression structure (can't compare
   cell-for-cell against the lat/lon result, since the grids differ).

## 8. Sequencing
1. **Finish Cartesian** — ENH-1 (flat resolution), ENH-2 (bowl-interior tiles), ENH-3 (internally
   drained), and the real-MPI driver. This maps the geometry interface.
2. **Geometry interface map** — enumerate every area/stencil/singularity touch-point (§5 discipline).
3. **Grid generation** — pick the construction (equal-area icosahedral leading), generate at target
   resolution.
4. **Port the geometry layer** — area accessor (leg a, mostly the WTM port), hex stencil + flowdir
   encoding (leg b), pentagon/seam handling (leg c), per-seam index transforms (§4).
5. **Validate** — analytic sphere surfaces + the serial differential oracle, exactly as for the Cartesian
   build.
