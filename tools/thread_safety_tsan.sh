#!/usr/bin/env bash
# thread_safety_tsan.sh -- run a ThreadSanitizer build of the distributed DepressionHierarchy shim and
# fail if TSan reports a real data race. This is the regression guard for concurrency bugs like the
# GDALAllRegister() heap-corruption race (rank threads constructing Array2Ds concurrently) that aborted
# ~1.5% of runs before it was fixed -- TSan detects that class DETERMINISTICALLY (on essentially every
# run), unlike a probabilistic stress loop.
#
# Two known non-bugs are removed AT THE SOURCE rather than via TSan suppressions (GCC libgomp's
# OpenMP-outlined functions don't match TSan's `called_from_lib`/`race:` suppression forms reliably, so a
# suppressions file here would be a false promise):
#   * OMP_NUM_THREADS=1 -- runs the flood's OpenMP regions serially, removing libgomp's spurious
#     "reduction race" reports (GCC libgomp isn't TSan-instrumented). This does NOT weaken the check we
#     care about: the real historical bug (the GDALAllRegister race) is in the RANK-thread concurrency,
#     fully exercised regardless of OpenMP. (Races *inside* an OpenMP region would need an
#     Archer-instrumented OpenMP runtime to test cleanly -- a known gap, not achievable with GCC.)
#   * RICHDEM_NO_PROGRESS (set on the build target) -- removes richdem's ProgressBar writing shared
#     std::cerr from every rank thread (a real, if benign, race).
# We run the flat-replay multi-tile case (most cross-rank machinery). The binary's exit code is ignored
# (it returns non-zero on a benign flowdir diff); pass/fail = "did the run finish (MPI-TREE-MATCH) with no
# ThreadSanitizer data race".
#
# Usage: thread_safety_tsan.sh <dephier_mpi_tsan.exe> <DEM> <ocean> <splits>
set -u
BIN="${1:?usage: thread_safety_tsan.sh <tsan-binary> <dem> <ocean> <splits>}"
DEM="$2"; OCEAN="$3"; SPLITS="$4"

# System libstdc++ first (dodge the toolchain/anaconda libstdc++ clash); keep the rest for libgdal etc.
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
export TSAN_OPTIONS="halt_on_error=1 detect_deadlocks=0"   # detect_deadlocks off: libgomp lock-order noise

# TSan and ASLR high-entropy layouts clash on some kernels ("unexpected memory mapping" fatal); run with
# ASLR disabled when setarch is available, else run directly.
runner=(); command -v setarch >/dev/null 2>&1 && runner=(setarch -R)

out="$( OMP_NUM_THREADS=1 DH_FLAT_PARTITION_REPLAY=1 "${runner[@]}" "$BIN" "$DEM" "$OCEAN" "$SPLITS" 2>&1 || true )"

if grep -q "ThreadSanitizer: data race" <<<"$out"; then
  echo "THREAD-SAFETY FAIL: ThreadSanitizer reported a data race"
  # Trim richdem progress bars, show the report.
  tr '\r' '\n' <<<"$out" | grep -vE '^\[|threads\)' | sed -n '/ThreadSanitizer: data race/,/^SUMMARY/p' | head -40
  exit 1
fi
if ! grep -q "MPI-TREE-MATCH " <<<"$out"; then
  echo "THREAD-SAFETY FAIL: run did not complete cleanly (no MPI-TREE-MATCH; crash or TSan fatal?)"
  tr '\r' '\n' <<<"$out" | grep -vE '^\[|threads\)' | tail -20
  exit 1
fi
echo "THREAD-SAFETY OK: no data race; run completed (MPI-TREE-MATCH)"
