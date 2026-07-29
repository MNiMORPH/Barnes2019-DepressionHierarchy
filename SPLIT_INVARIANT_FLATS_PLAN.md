# Plan: flag-gated split-invariant flats (GitHub #3 / ENH-7)

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
