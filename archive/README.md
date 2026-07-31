# Archived / historical documents

These are **superseded planning, design, engineering, and review records**, kept for provenance. They
describe earlier states of the project or completed investigations and are **not current documentation** —
so don't read them as guidance for how the code works today (they predate, e.g., the phase-function rename).

- `PARALLEL_DEPHIER_DESIGN.md` — early strategic scoping of the distributed build (dated, "for deciding
  direction"); the direction was chosen and built.
- `PARALLEL_DEPHIER_PLAN.md` — the plan the distributed build followed; now executed.
- `PARALLEL_DEPHIER_ENGINEERING.md` — the MPI-build engineering plan; now executed.
- `SPLIT_INVARIANT_FLATS_PLAN.md` — the split-invariant-flats investigation (self-marked SUPERSEDED); the
  flat-label replay it led to is ENH-8.
- `WTM_API_REVIEW.md` — a point-in-time external review of the DH public API, written before the
  phase-function rename (so it uses the old `PhaseAB/PhaseC/PhaseCD` names).

These docs use the **old phase-function names** (`GetDepressionHierarchyPhaseAB/PhaseC/PhaseCD`); the
authoritative old→new mapping is in `../DH_API_NAMING_REVIEW.md`.

**Current documentation** lives at the repository root: `../README.md`, `../VALIDATION.md` (how the
distributed build is validated), `../BENCHMARKS.md` (footprint + run-mode guidance), `../ENHANCEMENTS.md`
(the ENH log), and `../RICHARD_REVIEW_NOTES.md` (core changes for upstream review).
