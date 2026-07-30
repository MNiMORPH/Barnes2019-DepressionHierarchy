# DH API naming & granularity review — the fork's function surface, from the inside

Internal review of the function surface this fork grew while building the distributed DepressionHierarchy:
the **core phase functions** (`include/dephier/dephier.hpp`) and the **distributed-build tool helpers**
(`tools/dh_*.hpp`, `tools/comm_*.hpp`, and the two tool entry points). Companion to `WTM_API_REVIEW.md`,
which reviews the *core public API from the WTM consumer's outside view*; this one is the authors' inside
view, focused on **naming consistency and granularity** — the "functions have gotten a bit wild" concern.

Written by Claude with A. Wickert, 2026-07-30. Every line reference checked against source on 2026-07-30,
not recalled.

> **STATUS — rec 1 (core phase rename) DONE 2026-07-30** (commits `a890659` rename + `9ddd4c0` comment
> reframing). Authoritative old→new mapping (older docs still use the old names as historical record; this is
> the translation):
> `GetDepressionHierarchyPhaseAB` → **`FloodAndAssignDepressions`** (paper A+B+C, one grid pass);
> `GetDepressionHierarchyPhaseC` → **`ConstructHierarchy`** (paper D, grid-free primitive);
> `GetDepressionHierarchyPhaseCD` → **`ConstructHierarchyAndVolumes`** (central-only convenience);
> `GetDepressionHierarchy` unchanged. Paper letters now live in `dephier.hpp` comments where they apply.
> Recs 2–5 (tool-helper names, casing, comm symmetry, monolith extraction) not yet started.

## TL;DR

The **decomposition is right — the seams sit where they should**; it's the **names** that are wild. Four
things:

1. The core phase functions carry **shifted paper letters** — `PhaseAB` is actually the paper's A+B+**C**,
   `PhaseC` is the paper's **D**, `PhaseCD` tacks on volumes (§6.4, not a lettered phase). Rename to
   *behavior*, put the paper letters back in *comments*, and **do not re-carve** (see the granularity
   section — the fused flood is correct).
2. The shared tool headers **mix two casing conventions** (`CollapseSeamArtifacts`/`CommSend` vs
   `outlet_scan_intra`/`resolve_flat_flowdirs`/`flat_seed`).
3. A **six-member `resolve_flat_flowdirs*` family** whose suffixes encode *design history* ("option2")
   rather than the axis of variation. All six are live (verified); the issue is purely naming.
4. **Granularity of the tool entry points:** `rank_main` (506 lines, 37 inline lambdas) and
   `dephier_stitch`'s `main` (~605 lines) hide the very phase structure the core API exposes.

None of this is a correctness issue — it's legibility and future-proofing. Crucially, the core phase names
are **ours** (introduced this arc: `849578b`, `1714c53`; the 2020 original `b26a93b` had a single
`GetDepressionHierarchy`), so they are free to fix now and simply travel in the eventual upstream PR.

## Granularity: match the seams, not the paper's section headers

The paper (Barnes, Callaghan & Wickert 2020, §3) names four phases: A ocean identification, B pit-cell
identification, C depression assignment, D hierarchy construction (volumes are §6.4, post-processing). The
tempting move is one function per lettered phase. **That would be wrong.**

A function boundary earns its place when a **distinct intermediate object crosses it** and/or it is a
**parallel/serial or reuse boundary** — not when the paper has a heading there. Checked against source
(`dephier.hpp:285+`): `PhaseAB` gathers `ocean_seeds` (A) and `land_seeds` (B) into two vectors, then runs
**one** priority-flood loop (C) that mutates shared state — the PQ, the `outlet_database`, `label`,
`flowdirs`, `depressions` — together. A, B, C are three facets of one indivisible flood; you never call B
without A or run C on anything but the seeded PQ. Splitting them would be *false granularity*: three
functions threading the same mutable PQ, none independently callable, no intermediate product to show for
it.

Apply the test and there are exactly **two real seams → three stages** (which the code already has):

| Stage | Paper | Why it's a real boundary | Today | Proposed name |
|---|---|---|---|---|
| Flood + assign + outlet discovery | A+B+C | one grid pass; **per-tile parallel**; emits {leaf depressions, outlet set} | `GetDepressionHierarchyPhaseAB` | `FloodAndAssignDepressions` |
| Hierarchy construction | D | **grid-free, inherently central**; consumes outlets, emits the meta-tree | `GetDepressionHierarchyPhaseC` | `ConstructHierarchy` |
| Volumes | §6.4 | post-processing; marginal (distributable) + total (central) | `CalculateMarginalVolumes` / `CalculateTotalVolumes` | keep (already behavior-named) |

So: **rename, don't re-carve.** Put the paper's letters back in comments (`// Phase A (§3.1): seed ocean
cells` over the `ocean_seeds` block, etc.) so the code maps to the paper without pretending A/B/C are
separable functions. Leave A/B as inline seed-gathering (a dozen lines each, with OpenMP reduction pragmas —
two one-call helpers would read worse).

## The core phase functions (ours to rename; `GetDepressionHierarchy` is not)

| today | is really | rename to |
|---|---|---|
| `GetDepressionHierarchyPhaseAB` (`dephier.hpp:285`) | paper A+B+C (flood + outlet discovery) | `FloodAndAssignDepressions` |
| `GetDepressionHierarchyPhaseC` (`:653`) | paper D (hierarchy construction) | `ConstructHierarchy` |
| `GetDepressionHierarchyPhaseCD` (`:825`) | paper D + §6.4 volumes (central-only convenience) | `ConstructHierarchyAndVolumes` (keep) |
| `GetDepressionHierarchy` (`:849`) | Richard's one-call wrapper (preserved) | **leave alone** — his interface + the WTM drop-in target |

**Why the repeated stem is correct, not a smell (verified 2026-07-30).** `ConstructHierarchy` (`PhaseC`) is a
**primitive both paths use**; `ConstructHierarchyAndVolumes` (`PhaseCD`) is the **central path re-bundling**
that primitive with the two volume calls. The split is load-bearing: the distributed build (`dephier_mpi`)
calls `ConstructHierarchy` alone on rank 0 (`:792`), then *reimplements* the marginal-volume walk per-rank
and reduces it (`:804–818`) before a central `CalculateTotalVolumes` (`:819`) — it never calls the combo.
So `PhaseCD`'s only callers are the two **central/whole-grid** paths: the serial `GetDepressionHierarchy`
wrapper (`:857`) and the in-process `dephier_stitch` oracle (`:518`). The shared `ConstructHierarchy` stem
honestly names the shared primitive; keep the combo (it is not dead, and dropping it would only inline three
calls at two central sites for no gain). Its comment can say "centralized convenience" to signal it is not
the distributed path.

## The tool helpers

**The `resolve_flat_flowdirs*` family — six variants, all live, ad-hoc suffixes** (`tools/dh_flats.hpp`).
The suffix scheme mixes *method* ("option2" = the label-free relaxation), *distribution* (tiled / rank /
distributed), and *plumbing* (into), so a reader can't infer the axis of variation:

| fn | line | role | caller |
|---|---|---|---|
| `resolve_flat_flowdirs` | 33 | richdem-based reference (whole grid) | stitch serial ref (`:537`); test reference |
| `resolve_flat_flowdirs_option2` | 96 | label-free whole-grid relaxation (ENH-1) | `flat_mask_reconstruct_test` (equivalence witness) |
| `resolve_flat_flowdirs_option2_tiled` | 120 | per-tile relaxation + seam exchange (ENH-1 driver B) | `flat_mask_reconstruct_test` |
| `resolve_flat_flowdirs_rank` | 169 | per-rank (Comm) relaxation | `dephier_mpi` real flat pass (`:869`) |
| `resolve_flat_flowdirs_into` | 245 | sub-window helper | `_distributed` (`:271`) |
| `resolve_flat_flowdirs_distributed` | 262 | adaptive-halo capped (option 3) | stitch (`:317`) |

*Recommend:* drop the "option2"/"option3" history-jargon and name on a consistent axis — a base
`ResolveFlatFlowdirs` plus explicit distribution variants (e.g. `…Relaxed` for the label-free form,
`…Tiled`, `…PerRank`, `…AdaptiveHalo`), with `_into` demoted to a clearly-internal helper. (This is the
label-free ENH-1 sibling of the label-side `flat_label_distributed`.) All six are used, so this is renaming,
not removal.

**Casing is split across the tool headers.** PascalCase: `CollapseSeamArtifacts` (`dh_collapse.hpp:50`),
`CommSend`/`CommRecv`/`CommInit` (`comm_*.hpp`), the `OutletDB` type. snake_case: `outlet_scan_intra`/
`outlet_scan_seam` (`dh_outlets.hpp:57,77`), the whole `resolve_flat_flowdirs*` family, `flat_at`/`flat_seed`/
`flat_finish` (`dh_flats.hpp:52,57,72`), `canonicalize`/`invariants`/`relabel` (`dh_canonical.hpp`). The
**core library is PascalCase** (`GetDepressionHierarchy`, `CalculateMarginalVolumes`). *Recommend:* align
the DH-operation helpers to the core's PascalCase; the `dh_canonical` test-internal helpers are a judgment
call (test-local, lowercase is defensible). This one is genuinely a taste pick — flagging, not prescribing.

**The comm backends aren't a symmetric pair** (`comm_thread.hpp` / `comm_mpi.hpp`). Both expose
`CommSend`/`CommRecv`/`CommBarrier`/`CommRank`/`CommSize`, but init/finalize diverge: `comm_thread` has
`CommInit(n, fn)` (`:42`) only; `comm_mpi` has `CommInitMPI()` (`:33`) + `CommFinalizeMPI()` (`:34`) **and**
a stub `CommInit(n, fn)` (`:68`). So the driver's `#ifdef DH_USE_MPI` must know which init/finalize story
each backend tells. *Recommend:* a matched surface — the same `CommInit`/`CommFinalize` names in both, with
the MPI backend's extras (if any) documented as the deliberate difference.

## Granularity of the tool entry points (optional, lower priority)

`dephier_mpi.cpp`'s `rank_main` is **506 lines (374–880) with 37 inline lambdas** (`strip_of`,
`drain_local`, `localwalk`, `chase`, `pack`, `absorb`, `replay`, …); `dephier_stitch.cpp`'s `main` is
**~605 lines** (starts `:84`). The pipeline stages (per-rank flood, namespace remap, conduit resolution,
outlet set, flat-label reconciliation, tree gather, distributed volumes) are clearly *commented* but not
*named as functions* — the granularity lives in comments, not the type system, which is ironic given the
core-API argument above. *Recommend (only if/when it's worth the churn — correctness is done):* extract the
named stages into functions, mirroring the very phase decomposition the core exposes. This is readability,
not behavior; do it deliberately, not as a drive-by.

## Recommendations, ranked

1. **Rename the core phase functions to behavior; paper letters to comments; do not re-carve.**
   `PhaseAB`→`FloodAndAssignDepressions`, `PhaseC`→`ConstructHierarchy`, `PhaseCD`→
   `ConstructHierarchyAndVolumes` (keep — the central-only convenience over the shared `ConstructHierarchy`
   primitive; see the phase table's note). Keep `GetDepressionHierarchy`, keep the volume names. These
   names are ours (`849578b`/`1714c53`), so free to do now; they ride the eventual Richard PR for his
   blessing, with no "rename his interface" hazard.
2. **Rename the `resolve_flat_flowdirs*` family on a consistent axis; drop "option2/3".** All six live;
   naming only.
3. **Pick one casing convention** (recommend the core's PascalCase for DH-operation helpers).
4. **Symmetrize the comm backends** (matched `CommInit`/`CommFinalize`).
5. **(Optional, lower)** Extract the tool-`main` pipeline stages into named functions.

## Verification (checked against source 2026-07-30)

- Core phases: `GetDepressionHierarchyPhaseAB` `dephier.hpp:285` (gathers `ocean_seeds`/`land_seeds` then one
  flood loop — A/B/C fused over shared state), `PhaseC` `:653`, `PhaseCD` `:825`, `GetDepressionHierarchy`
  `:849`. Phase functions introduced by `849578b` (AB/CD split) and `1714c53` (C/CD split); original
  `b26a93b` (R. Barnes, 2020-04-16) had a single `GetDepressionHierarchy`, no phase functions, no A/B/C/D.
- Flat family: `dh_flats.hpp:33,96,120,169,245,262`; helpers `flat_at:52`, `flat_seed:57`, `flat_finish:72`,
  `NOT_FLAT:244`. Call sites: `dephier_stitch.cpp:537,317`, `dephier_mpi.cpp:869`,
  `flat_mask_reconstruct_test.cpp:111,112,114`, `_distributed`→`_into` at `dh_flats.hpp:271`.
- Casing: `CollapseSeamArtifacts` `dh_collapse.hpp:50`; `outlet_scan_intra/seam` `dh_outlets.hpp:57,77`;
  `CommSend/Recv/Init/Barrier/Rank/Size` `comm_thread.hpp:52,61,42,79,37,38`.
- Comm backends: `comm_thread.hpp` `CommInit:42` (no finalize); `comm_mpi.hpp` `CommInitMPI:33`,
  `CommFinalizeMPI:34`, stub `CommInit:68`.
- Monoliths: `dephier_mpi.cpp` `rank_main` 374–880 (506 lines), 37 inline `const auto … = [` lambdas;
  `dephier_stitch.cpp` `main` from `:84`, file 689 lines.
