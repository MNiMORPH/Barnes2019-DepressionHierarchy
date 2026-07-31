# Validating the distributed DepressionHierarchy

How do we know the tiled / distributed build is correct — including on inputs (GEBCO 30″) far too large to
check against a serial run? The answer is two complementary pillars: **exactness against serial where serial
can run, and split-invariance where it cannot.** The first proves the algorithm is right; the second proves
it stays right at scales the first can't reach.

## Pillar 1 — Exactness: an equivalence proof against serial

The distributed build can reproduce the **bit-identical** serial depression hierarchy — same tree, down to
the canonical signature (per-node `pit_elev`, `out_elev`, `cell_count`, `dep_vol`; see `tools/dh_canonical.hpp`).
This is not a spot-check; it is an equivalence demonstration of the *entire* pipeline:

- per-tile flood (`FloodAndAssignDepressions`),
- namespace remap (per-tile → global depression labels),
- conduit resolution (BOUNDARY cells resolved across seams),
- the flat-label replay (`DH_FLAT_PARTITION_REPLAY`, uncapped) — reproduces serial's exact flat partition,
- outlet re-derivation + cross-seam `HandleEdge`,
- grid-free hierarchy assembly (`ConstructHierarchy`),
- distributed volumes (per-rank marginal walk reduced to rank 0).

If all of that yields serial's exact tree, the decomposition has not changed the computation — the network
*is* serial's algorithm, just cut into tiles and message-passed. Verified across the committed fixture suite
plus exhaustive single-seam and multi-seam fractal sweeps (adversarial fractals, Corsica, pit-on-seam
fixtures), encoded as bit-identical CTest cases: `mpi_flat_partition_replay`, `mpi_flat_replay_chained`,
`mpi_flat_replay_canonical`, `flat_partition_replay_proof`, and the `stitch_*` in-process checks. The
in-process stitch matches serial bit-for-bit across the whole sweep; the real-MPI build validates the tree
on rank 0 (`mpi_real_tree_*`).

Bit-identity rests on three deterministic-flood fixes that **change serial output** and are therefore flagged
for upstream review (see `RICHARD_REVIEW_NOTES.md`): the `ConstructHierarchy` geometric outlet tie-break, the
`radix_heap` equal-elevation ordering, and the `out_cell` tie-break. It also relies on *reproducing* serial's
flat labelling (the pit-index flood replay) rather than altering it — a match to the algorithm, not a change.

### What this buys for the non-exact run modes

The exact mode is opt-in. The default (no replay) and capped-replay modes are **not** unvalidated code — they
are the *same proven pipeline* run with a deliberate, characterized relaxation. Every difference from serial
falls in a named, measured class — flat-cell tie-breaks and the meta-vs-`ocean_linked` binarization tie-break
— and every one is `VOL-MATCH` (volume-correct) and a valid tree. So the non-exact modes are a validated
engine at a bounded, understood approximation, not a leap of faith. (See `BENCHMARKS.md` → "Choosing a mode".)

## Pillar 2 — Split-invariance: a serial-free self-check at scale

At GEBCO 30″ (~272 M land cells) serial cannot run — it does not fit in one rank and is far too slow — so the
Pillar-1 comparison is unavailable. Confidence there rests on:

1. **The machinery is proven equivalent on representative, tractable inputs**, and the design is
   size-agnostic (no step's correctness depends on the grid fitting in memory; footprint is `O(N/P) +
   O(boundary)` — see `BENCHMARKS.md`). Correctness demonstrated at testable scale carries to larger N.
2. **Split-invariance** — the direct self-check that needs no serial reference. A correct build produces the
   **identical tree regardless of where the seams fall**, so:
   - `STITCH-SIG` (a hash of the build's own canonical tree) must be identical across two different tilings of
     the same DEM — "same DH no matter how we split," which references only the build itself; and
   - `*-DECOMP-CORRECT` (node count independent of the split) is a necessary structural invariant — an unequal
     node count between tilings is a definitive decomposition error.

   These catch decomposition bugs (a seam that changed the tree) without ever needing to run serial.

## Summary

- **Exactness vs. serial** — where serial runs — proves the distributed pipeline computes serial's tree
  exactly: an equivalence proof, not a plausibility check.
- **Split-invariance** — where serial cannot run — proves the build stays internally consistent (same tree at
  any tiling) at the scales that motivate distributing in the first place.

Together they are why we can trust the network on inputs serial will never see: the algorithm is proven
correct where it's checkable, and proven stable where it isn't.
