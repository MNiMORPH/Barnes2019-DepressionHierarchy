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

## Approach
In flag mode, compute a **canonical floor-flat labelling on the full grid** (geometry only, tiling-
independent — acceptable under the reproducibility/speed trade; later distributable via the same ENH-1
seam-exchange relaxation). A connected floor-flat (is_flat = no strictly-lower D8 neighbour; sills excluded)
is ONE depression whose canonical pit = its **lowest-global-index cell**. Merge every tile-half leaf whose
pit lies in that flat into one canonical leaf, relabel all its cells, compact away the emptied leaves, then
run the existing outlet re-derivation + PhaseCD + collapse on the canonicalised labels/tree.

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
