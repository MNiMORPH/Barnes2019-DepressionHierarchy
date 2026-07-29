# Plan: flag-gated split-invariant flats (GitHub #3 / ENH-7)

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
