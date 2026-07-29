# Plan: flag-gated split-invariant flats (GitHub #3 / ENH-7)

> **⚠ SUPERSEDED / HISTORICAL (2026-07-30).** This is the working record of the split-invariant-flats
> investigation, including two approaches that were tried and are now REMOVED from the code: *source-1*
> (`DH_SPLIT_INVARIANT_FLATS`, subsumed by the replay) and the *containment cc-pass*
> (`DH_CONTAINMENT_CCPASS`, validated-negative — it computes drainage, but serial's flat label is flood
> membership). The approach that shipped is the flat-partition **replay** (`DH_FLAT_PARTITION_REPLAY`), now
> the primary path, distributed per-rank in `dephier_mpi` as **ENH-8**. Current state lives in
> `ENHANCEMENTS.md` (ENH-7/ENH-8) and `RICHARD_REVIEW_NOTES.md`; kept here as the dead-end record.

> **FULL IDENTITY IS ONE OUTLET-TIE-BREAK AWAY — root cause found + PROVEN (2026-07-29).** After the replay
> closed the cell-assignment class (73→79/107 MATCH), characterized the 28 residual DIFFERs: all are
> STITCH-LEAFSET-MATCH (committed `451073d`) — exactly serial's leaf depressions, difference is pure meta-tree
> SHAPE (8 meta-ordering + 20 meta-vs-ocean_linked). Root cause (DH_AUDIT_VS_SERIAL, `d80c41d`): PhaseC's sort
> is deterministic on (out_elev → **out_cell** → pit cells); the tiled outlets match serial in pairs + out_elev
> but differ in out_cell — serial keeps the FIRST-DISCOVERED out_cell on a tie (flood-order, dephier.hpp:579),
> the tiled scan keeps the LOWEST (tiling-independent, dh_outlets.hpp). At a tied out_elev this flips the sort.
> **PROVEN (experiment, reverted): replay + a 1-line serial change to prefer the lowest out_cell on a tie =>
> 107/107 STITCH-MATCH.** The change is serial-output-changing (Richard) but makes serial itself
> split-invariant, and Barnes' own comment (dephier.hpp:42) already calls the tie "arbitrarily chosen." This is
> lever B, now precisely scoped. Lever C = distribute the replay (ENH-1 seam-exchange). Old "lower cell wins
> reproduces serial" claim is DISPROVEN — serial keeps flood-first, not lowest.

> **ROAD A — flat-label PARTITION replay: rule DERIVED from the flood + (partially) VALIDATED (2026-07-29).**
> Andy chose Road A (reproduce serial's flat-label partition locally, serial untouched). Derived the exact
> rule from the flood (`dephier.hpp:481-536`) + the fork's radix_heap tie-break (`radix_heap.hpp:391-399`:
> bucket-0 sorted ASC by cell index, consumed `pop_back` ⇒ equal-elevation cells pop HIGHEST-INDEX first;
> same-elevation cells labelled mid-bucket are appended and popped LIFO). Rule for a connected same-elevation
> region: every lower bucket drains before bucket e, so each region cell with a strictly-lower neighbour is
> SEEDED (before bucket e) by the label of its lowest-(elev), highest-(index) strictly-lower neighbour —
> **OCEAN counted at its dem elevation** (an ocean cell is a seed popped at its dem value; this is the ocean
> handling the earlier sill attempt got wrong); bucket e then propagates into the `is_flat` interior
> (no strictly-lower neighbour) HIGHEST-INDEX-first + LIFO. A closed floor-flat (no lower exit) = a pit,
> highest-index cell seeds a fresh basin (= source-1's domain).
>
> **CLEAN PROOF DONE (`scratchpad/flat_full_replay.cpp`): the rule reproduces serial's label partition
> EXACTLY over EVERY land cell — 49 DEMs, 0 partition-differ.** (24 kerry/testdem fixtures + Corsica 14 064
> land cells + 24 adversarial fractals, sizes 120/200, β 1.3–2.2, seeds 1/4/7/11.) The proof is a full
> pit-up replay faithful to the flood: seeds = every OCEAN cell + every land cell with NO strictly-lower
> neighbour (the land_seed / pit set, `dephier.hpp:412` — INCLUDES all is_flat interiors), each pushed at
> its dem elevation; process buckets elev-ASC, within a bucket sort ASC by index and consume `pop_back`
> (highest-index first) with same-elevation labels appended + popped LIFO; a popped NO_DEP cell becomes a
> new pit. Compared namespace-free (bijection s_label↔replay) over all land cells. The KEY correction over
> the first isolated test: ALL is_flat cells are land_seeds already in the PQ, so a flat interior pops as a
> NEW pit only if still NO_DEP when popped (highest-index-first) — the first test's "seed-then-propagate"
> missed this, producing the closed-floor-flat confound; the full replay eliminates it. So the flat-partition
> rule + dynamics are now FULLY understood and TESTED, with a correct reference implementation to build from.
>
> **INTEGRATION (i) BUILT + VALIDATED (2026-07-29): `DH_FLAT_PARTITION_REPLAY`.** Full-grid replay (§above)
> overwrites gLabel with serial's partition mapped onto G's existing leaves (each replay-basin → the min-label
> G-leaf whose pit falls in it; several leaves in one basin = a seam-split → they merge, SUBSUMING source-1),
> stamping each canonical leaf with the replay's OWN pit cell (= serial's pit, so pit_elev — which IS in the
> signature — matches; fixed testdem8 sp3). Then COLLAPSE is RETIRED under the flag: the correct partition
> means outlets→PhaseCD build serial's tree with no artifacts and the compact drops emptied leaves, so the
> collapse's seam-dependent meta-over-halves would only DISSOLVE REAL structure (measured: kerry_test9 sp7
> merged two genuine basins until the collapse was skipped). **Sweep (107 fixture×split cases): STITCH-MATCH
> 73→79 (+6, ZERO broke), 0 VOL-DIFFER, 0 orphan cells, 0 DECOMP regressions; default OFF byte-identical
> (suite 23/23).** The residual DIFFERs are now the SEPARATE levers: pure meta-nesting/outlet-ordering (DIFFER
> but DECOMP-CORRECT — the leaf cell_counts already match serial, e.g. kerry_test5 sp3) and meta-vs-ocean_linked
> (DECOMP-INCORRECT) — both Richard-coordinated, NOT cell-assignment. The cell-assignment DIFFER class is
> CLOSED by the partition fix. Open design Qs (for Andy): (1) the replay subsumes source-1 — unify/deprecate
> the DH_SPLIT_INVARIANT_FLATS flag? (2) collapse retirement is now real under replay — path to default? (3)
> the superseded containment cc-pass code (DH_CONTAINMENT_CCPASS + collapse out_map) — remove now?
>
> **INTEGRATION DECISION (proof now DE-RISKS option i):** (i) run the replay full-grid over gLabel
> (source-1 style, flag-gated) — now JUSTIFIED: proven to reproduce serial's partition exactly, so it cannot
> break non-seam cells; fastest way to MEASURE whether the partition fix (with source-1) collapses the
> cell-assignment DIFFER class in the stitch. Distributable later, same status as source-1. (ii) target ONLY
> seam-straddling flats via ENH-1 seam-exchange — the honest distributed form; the propagation is sequential
> (highest-index LIFO, new-pit-on-pop), so it replays the flood ORDER across the seam (bounded per flat,
> global-row-major index is consistent across tiles) rather than independent relaxation rounds like ENH-1's
> flat FLOWDIRS. RECOMMEND: do (i) in the stitch next to confirm the identity gain, then (ii) as the
> distributed engineering. Tools kept: `flat_full_replay.cpp` (the proof/reference), `flat_replay_test.cpp`.

> **FINDING (2026-07-29) — the containment cc-pass via flowdir-tracing is the WRONG mechanism; the residual
> is a FLOOD-MEMBERSHIP tie on flats.** Built the cc-pass (flag `DH_CONTAINMENT_CCPASS`, post-collapse,
> structure-preserving: re-derive each cell's leaf by tracing the serial-identical fixed flowdirs `gFix` down
> to a pit, then recompute ONLY marginal+total volumes on the FROZEN tree; collapse extended with an optional
> pre→post label `out_map` to carry dissolved-pit cells to their surviving leaf). It did NOT reconcile:
> kerry_test9 sp7 meta 35→36 (serial 29); kerry_test10 sp5 fixed the top meta (195→**197** ✓) but **EMPTIED a
> zero-height flat leaf** `(6,6,20,0)→(6,6,0,0)` and mis-moved others `(8,8,49)→(8,8,26)`. VOL stayed MATCH
> (structure-preservation held — no NaN), but SIG did not converge.
>
> **ROOT (mechanistic, decisive):** a cell has TWO different "assignments" that COINCIDE on strict slopes but
> DIVERGE on flats: (1) DRAINAGE — where its flowdir sends water (a sill/through-flat's flat-flowdirs point at
> the EXIT, out to a lower basin); (2) MEMBERSHIP — which depression Priority-Flood LABELS it into (the basin
> whose wavefront claims it first, by (elevation, index) pop-order). Serial's leaf label = MEMBERSHIP.
> Flowdir-tracing computes DRAINAGE. So the cc-pass empties every zero-height flat leaf `(e,e,·,0)`: its cells
> drain out (→ neighbour pit) but serial COUNTS them in the flat. This is the SAME rock the source-2 per-cell
> steepest-descent hit ("split sill-flats, kerry_test2 cc 28/8 vs 32/4"), now proven from the attribution
> side. And attribution can't sidestep it: for a sill cell `my_elev==out_elev`, the walk-up `while(my_elev >
> out_elev)` does NOT fire, so the cell is pinned to whatever leaf it's LABELLED with — the partition directly
> sets the zero-height leaves' cell_count. No attribution trick recovers it; you must match serial's flat
> partition. **Kept as a measured negative; not reverted yet (Andy).** Diagnostics added and worth keeping:
> `DH_DUMP_FDDIFF` (flowdir divergences — 100% seam-confined, mostly flats) and `DH_CCPASS_CEIL` (gFix-trace
> vs serial-trace terminal-pit mismatch = the drainage residual; 0 on kerry_test9 sp7, i.e. drainage was
> perfect there yet SIG still differed → proof the block is membership/structure, not drainage).
>
> **DEEPER APPROACH TO RECONCILE (the real lever):** reproduce serial's deterministic flat-LABEL PARTITION,
> not a flowdir trace. Source-1 already handles the trivial case (a FLOOR flat — all cells below their
> neighbours, no competition → one basin). The unsolved case is the SILL/through-flat: cells at their own
> elevation with equal-elevation exits to ≥2 basins, which serial partitions by flood pop-order (a geodesic
> Voronoi from the competing spill cells, tie-broken by index). Roads: **(A)** compute that partition locally
> per connected flat — deterministic multi-source BFS from the flat's basin-adjacent perimeter, ordered by
> (distance, global cell index) to match the flood, distributable via ENH-1 seam-exchange on the LABEL field
> (this is "source 2 as a label relaxation," the principled form; the earlier sill attempt failed only
> because its rule — lowest-adjacent-lower, ocean as −∞ — was the WRONG partition, not because partitioning is
> wrong). **(B)** canonicalize the tie in BOTH serial and tiled (Richard co-scope): change the DH core's flat
> labelling to a deterministic rule (e.g. lowest-index basin wins) so tiled matches by construction. (A) keeps
> serial untouched but must replicate its exact pop-order; (B) is a coordinated upstream change that makes the
> tie well-defined. Either way the mechanism is LABEL-PARTITION, and the collapse `out_map` hook stays useful.

> **GOAL RAISED (Andy, 2026-07-29):** aim for FULL serial↔parallel IDENTITY (bit-identical canonical
> signature), not merely volume-correct or split-invariant. PRINCIPLE (verbatim, "definitely remember"):
> **"They're complementary and both necessary: the collapse gets the right set of depressions; source 2
> gets the right cells into each."** For IDENTITY, each tie-break must reproduce SERIAL's specific choice.
>
> **COMMITTED & SAFE:** source-1 early floor-flat unification (`010c4a4`, flag DH_SPLIT_INVARIANT_FLATS) +
> collapse meta-over-halves retirement warning (`3e28f7c`). Default OFF byte-identical (suite 23/23). The
> flag drives the seam-dependent meta-over-halves collapse (Pass B/B2) to 0 on floor-straddle flats. Also
> committed this arc: STITCH-SIG split-invariance metric + `tools/check_split_invariance.sh` (`98c7613`).
>
> **FULL SPLIT-INVARIANCE MAP (flag OFF):** 8/28 fixtures SPLIT-VARIANT, ALL volume-correct:
> kerry_test 1,2,4,5,9,10,11,12. Source-1 flag alone makes NONE of them invariant (variance dominated by
> the cell-assignment / structure levers below, not the collapse-work retirement).
>
> **SOURCE-2 ATTEMPTS — BOTH BROKE BROAD (the recurring trap):**
> - per-cell steepest-descent (scratchpad `dephier_stitch_source2_percell_wip.cpp`): split sill-flats →
>   kerry_test2 split-invariant but SERIAL-WRONG (cc 28/8 vs 32/4).
> - sill-region = lowest-adjacent-lower (scratchpad `dephier_stitch_source1_plus_sill_wip.cpp`): kerry_test2
>   RIGHT (cc 32/4, split-invariant) BUT (a)-sweep broke 11 bit-identical cases + corrupted volume in 11
>   (testdem8 → NaN). ROOT: treated NoData-ocean as −∞; a basin RIM adjacent to NoData-ocean got sent to
>   ocean, severing the basin → open depression → NaN. Serial keeps NoData-ocean at its SENTINEL elevation.
>
> **GENERAL SOLUTION (designed, Andy-informed, NOT yet built):** the relabel-before-outlets lets a
> labelling error ripple into STRUCTURE (a mis-labelled cell changes outlets → tree → NaN). So SEPARATE:
> (1) keep the collapse-corrected tree as the STRUCTURE (do not re-derive outlets from relabelled cells);
> (2) fix CELLS as a STRUCTURE-PRESERVING pass — recompute only marginal volumes (cell_count/dep_vol) on the
> frozen tree, assigning each cell to its CONTAINING leaf = the leaf it reaches by DETERMINISTIC drainage
> (ENH-1 flat flowdirs + steepest descent), with OCEAN AT ITS FLOOD ELEVATION (not −∞). Frozen tree ⇒ a
> labelling error can only mis-COUNT, never spawn a NaN basin. One containment rule subsumes floor/sill/rim;
> reuses ENH-1 flowdirs. LEAVES the meta-vs-ocean_linked difference (kerry_test2 residual) as a SEPARATE
> third lever (a PhaseCD binarization tie), untouched by the cell pass.
>
> **VERIFICATION PROTOCOL (Andy's "a then b", earned the hard way):** before trusting ANY source-2 attempt,
> run the (a) broad sweep first — flag ON vs true-OFF (`env -u`!), require broke-a-match=0, VOL-DIFFER=0,
> 0 crashes across all fixtures×splits. Only then chase serial identity (b). Narrow-green (kerry_test2) ≠
> done. Every geometric shortcut so far passed narrow and broke broad — get the sweep before claiming.
>
> **KEY INSIGHT — VOLUME IS COMPUTED FROM STRUCTURE (the reconciliation, the strongest argument):**
> `dep_vol = cell_count·out_elev − total_elevation`. So a STRUCTURAL break surfaces AS a volume anomaly: the
> testdem8 NaN was an empty OPEN depression `(4,inf,0,-nan)` — my rim→ocean relabel SEVERED the basin's spill
> path, orphaning a leaf with no outlet (out_elev=inf) and 0 cells → `0·inf = NaN`. Consequences that drive
> the design: (i) FREEZING the tree makes the NaN class IMPOSSIBLE (every out_elev finite ⇒ dep_vol finite);
> (ii) `VOL-MATCH` is a STRUCTURE check, not just a volume check — it caught the break instantly, so it is the
> load-bearing (a)-sweep tripwire; (iii) with a frozen tree a labelling error can only shift a cell's marginal
> volume between finite-out_elev leaves → a *finite* VOL-DIFFER (caught, recoverable), NEVER a NaN. That
> downgrade — "broken tree" → "caught miscount" — is the entire reason source-2 must NOT touch structure.
> CONTAINMENT rule for the cc-pass: a cell → the leaf it reaches by DETERMINISTIC drainage (ENH-1 flat
> flowdirs + steepest descent), OCEAN AT ITS FLOOD/SENTINEL ELEVATION, not −∞ (the −∞ shortcut is exactly
> what sent basin rims to ocean).
>
> **PROCESS NOTE:** a premature "success ✓" is a destabilising false anchor; state what was tested AND what
> wasn't; hold green PROVISIONAL until verified at the real scope; ask "what would break this?" and check it
> FIRST so reversals die in my verification, not Andy's model. And: a partial step is a KEPT rung — don't
> reset to zero when it's not the whole staircase.


## Goal
An **opt-in** mode (flag) in which the tiled DepressionHierarchy is **split-invariant**: `STITCH-SIG`
identical across tilings for every DEM. Trades speed + complexity for reproducibility/safety (Andy's call).
Acceptance:
- kerry_test2/10/11 (currently SPLIT-VARIANT) become SPLIT-INVARIANT (`tools/check_split_invariance.sh`).
- kerry_test3/7/12 stay SPLIT-INVARIANT; all current `STITCH-MATCH` cases stay MATCH.
- **Flag OFF ⇒ byte-identical to today** (suite 23/23 unchanged).

## Root cause (established, two failed experiments)
A flat depression **floor** cut by a tile seam becomes a **multi-pit** floor-flat: each tile-half is a
separate leaf with its own flood-chosen (seam-dependent) pit cell. Downstream the two halves are unified by
outlets + PhaseCD + collapse, but the unified result (and a boundary cell shared with an adjacent
depression) stays seam-dependent → the tree varies with the seam. A label post-pass can't fix it: it must
give the flat ONE seam-independent identity and MERGE the split leaves.

## Approach — REDIRECTED to ENH-1 (2026-07-29, after the post-hoc surgery hit whack-a-mole)
The post-hoc full-grid gLabel surgery (S1-S4, saved `scratchpad/dephier_stitch_S1-S4_wip.cpp`) fixed
kerry_test2/3/7 but **regressed kerry_test12 and left 10/11** — late label-rewriting collides with the
outlet stage. Redirect to ENH-1's proven seam-exchange relaxation, done EARLY (before PhaseCD/collapse).

**Why early / why ENH-1 (grounded in dephier.hpp):** the cell label grid holds ONLY leaf labels; metas are
tree nodes, and `CalculateMarginalVolumes` (dephier.hpp:878) attributes each cell by walking UP from its
leaf label on the fly — the grid is never written with a meta label. So a label-carrying relaxation stays
in the leaf namespace and never meets a metadepression label. "Abandoned" labels are only real AFTER the
collapse dissolves/merges nodes (pre-collapse gLabel then points at stale labels) — running the relaxation
EARLY (live leaf namespace, before PhaseCD/collapse) sidesteps that entirely. Bonus: a cell at meta level
(≥ both children's outlets) is attributed to the meta by the walk-up regardless of leaf label, so the
relaxation need only nail BELOW-outlet leaf cells, and a closed floor-flat belongs to exactly one leaf.

**Two mechanisms (both via ENH-1's 1-column seam-exchange, so distributed + split-invariant by construction):**
1. **Floor-flat label unification** — a min-relaxation of the leaf label over same-elevation floor-flat
   adjacency: every connected floor-flat converges to ONE (min) label, seam-independently. This is S2 done
   the ENH-1 way (no full-grid flood-fill, no seam-dependent canonical choice). Merge the tile-half leaves
   accordingly; compact.
2. **Deterministic resolution for BOUNDARY / slope cells** — the conduit resolution currently follows the
   FLOOD's flat flowdirs (`tile.fd`, seam-dependent). Point it at ENH-1's DETERMINISTIC flat flowdirs
   (`resolve_flat_flowdirs`) instead, so BOUNDARY + slope cells resolve seam-independently. This replaces
   the failed post-hoc S4 descent trace with a principled "use the deterministic flowdirs we already trust."

Scope the exact touch-points (does the conduit read `tile.fd` where we can swap in the resolved flat fd?)
before coding; the last approach failed for lack of that grounding.

## Incremental steps (each measured with check_split_invariance.sh; flag OFF must stay 23/23)
- **S1 — region finder (build + eyeball).** Full-grid connected floor-flats + canonical (lowest-index) pit.
  Report region count / sizes on the 3 fixtures. No behaviour change yet.
- **S2 — merge + compact.** Build old-leaf→canonical-leaf remap for every leaf whose pit is in a region;
  set the canonical leaf's pit = the region's lowest-index cell; relabel gLabel; compact G to drop the
  merged-away leaves. Apply BEFORE outlet re-derivation.
- **S3 — measure.** STITCH-SIG across splits 2/3/5/7 on kerry_test2/10/11. Expect ≥ partial collapse of the
  distinct-tree count. Record the number; if variance remains, localise the residual cell(s).
- **S4 — residual.** If a boundary cell with an adjacent depression still flips, handle it (likely the
  outlet-cell tie on the flat rim resolved against the canonical pit). Iterate to SPLIT-INVARIANT.
- **S5 — productionise.** env flag → real CLI flag (e.g. `--split-invariant-flats`); confirm OFF is
  byte-identical (suite 23/23); add `split_invariance_kerry2/10/11` CTests asserting invariance ON; commit
  granularly. Update ENHANCEMENTS.md ENH-7 + issue #3.

## Guardrails
- Flag OFF is the default and must not change any output.
- Verify with the serial-free metric (check_split_invariance.sh), NOT stitch-vs-serial (a valid tie-break
  vs serial is fine; only cross-tiling consistency matters).
- Commit only after S5 validates; keep experiments env-gated and reverted until then.
