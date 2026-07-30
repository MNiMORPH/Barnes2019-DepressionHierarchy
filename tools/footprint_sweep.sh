#!/usr/bin/env bash
# footprint_sweep.sh -- quantify the per-rank footprint of the distributed DepressionHierarchy
# flat-label replay (DH_FLAT_PARTITION_REPLAY) and its footprint/exactness trade-off vs the halo cap.
#
# WHY: the flat replay's ADAPTIVE (uncapped) halo grows until each owned cell's replay label is
# stable, i.e. until the window reaches that cell's basin pit -- which on real terrain can be most of
# the grid. So the uncapped mode is bit-identical to serial but NOT footprint-bounded. A finite halo
# cap (4th CLI arg) bounds the footprint to owned_cols + 2*cap, but then only reproduces a VALID,
# volume-correct tree -- and bit-identity to serial is NON-MONOTONIC in the cap (a bigger halo is not
# reliably "more serial-identical"; see Sweep C). Sweeps: A = footprint/exactness vs cap; B = footprint
# vs tile count; C = exactness vs cap across terrain roughness (the non-monotonicity).
#
# Metric: DH_HALO_DIAG prints, per rank, owned_cols (= W/ntiles, the O(N/P) baseline) and held_cols
# (owned + halo = the peak transient footprint of the replay). We report the PEAK held_cols across
# ranks, plus whether the tree stayed bit-identical (MPI-TREE-MATCH) and the flowdir diff.
#
# Requires: build/dephier_mpi.exe (thread-shim build). The DEM defaults to a freshly generated 256^2
# fractal via make_synthetic_dem.py (needs WTM_DIR; see that script) -- or pass your own DEM as $1.
#
# Usage:   tools/footprint_sweep.sh [DEM] [BINARY]
# Example: tools/footprint_sweep.sh                 # generate + sweep a 256^2 fractal
#          tools/footprint_sweep.sh my.dem          # sweep an existing DEM
#
# (No `set -e`: the parsing greps legitimately return non-zero when a field is absent, and a missing
# value is reported as "?" rather than aborting the sweep.)

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${2:-$REPO/build/dephier_mpi.exe}"
# Put the system libstdc++ first (dodges the anaconda/toolchain libstdc++ clash); keep any existing
# path after it so a GDAL from e.g. anaconda is still found.
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"

DEM="${1:-}"
if [[ -z "$DEM" ]]; then
  DEM="$(mktemp --suffix=.dem)"
  echo "# generating a 256x256 beta=1.5 fractal DEM -> $DEM"
  python3 "$REPO/tools/make_synthetic_dem.py" --size 256 --beta 1.5 --seed 1 -o "$DEM" >/dev/null
fi
[[ -x "$BIN" ]] || { echo "error: $BIN not built"; exit 1; }

W=$(awk 'NR==1{print $2}' "$DEM")
echo "# DEM=$DEM  W=$W  binary=$BIN"

# Evenly-spaced interior split columns for N tiles, e.g. splits 4 -> "64,128,192" (W=256).
splits(){ local n=$1 i; local out=""; for ((i=1;i<n;i++)); do out+="$(( W*i/n )),"; done; echo "${out%,}"; }

# Run one config; echo "held_peak owned_max tree fd".
run(){ # args: split_string  cap(optional)
  local spec="$1" cap="${2:-}" out
  out=$(DH_FLAT_PARTITION_REPLAY=1 DH_HALO_DIAG=1 "$BIN" "$DEM" 0 "$spec" $cap 2>/dev/null || true)
  local held owned tree fd
  held=$(grep -oE 'held_cols=[0-9]+' <<<"$out" | grep -oE '[0-9]+' | sort -n | tail -1)
  owned=$(grep -oE 'owned_cols=[0-9]+' <<<"$out" | grep -oE '[0-9]+' | sort -n | tail -1)
  tree=$(grep -oE 'MPI-TREE-(MATCH|DIFFER)' <<<"$out" | head -1)
  fd=$(grep -oE 'fd_diff=[0-9]+' <<<"$out" | head -1)
  echo "${held:-?} ${owned:-?} ${tree:-?} ${fd:-?}"
}

NT=8; SPEC=$(splits $NT)
echo
echo "## Sweep A -- halo cap at fixed ntiles=$NT (owned_cols=$((W/NT)))"
printf "%-10s %-16s %-16s %-16s %s\n" "cap" "held_cols(peak)" "held/W" "tree" "flowdir"
for cap in "" 64 32 16 8 4; do
  read -r held owned tree fd < <(run "$SPEC" "$cap")
  printf "%-10s %-16s %-16s %-16s %s\n" "${cap:-inf}" "$held" "$(awk "BEGIN{printf \"%.2f\",$held/$W}")" "$tree" "$fd"
done

echo
echo "## Sweep B -- ntiles at fixed halo cap=16 (footprint = owned_cols + 2*cap)"
printf "%-10s %-16s %-16s %-16s %s\n" "ntiles" "owned_cols" "held_cols(peak)" "tree" "flowdir"
for nt in 2 4 8; do
  read -r held owned tree fd < <(run "$(splits $nt)" 16)
  printf "%-10s %-16s %-16s %-16s %s\n" "$nt" "$owned" "$held" "$tree" "$fd"
done

# Sweep C -- how bit-identity (MPI-TREE-MATCH) vs the cap varies with terrain roughness (beta). Generates
# its own beta-varied fractals, so it needs make_synthetic_dem.py (WTM_DIR); skipped if unavailable.
# Legend: M = MATCH (bit-identical); D = DIFFER but VOL-MATCH (valid, volume-correct, not serial-identical);
# ! = DIFFER *and* VOL-DIFFER (a real problem -- should never appear).
echo
echo "## Sweep C -- exactness vs cap across terrain roughness beta, ntiles=8 (M=match, D=differ-but-valid, !=bad)"
if python3 "$REPO/tools/make_synthetic_dem.py" --help >/dev/null 2>&1; then
  SPEC8=$(splits 8)
  printf "%-16s %s\n" "beta" "cap:  inf  64  32  16   8   4   2"
  for b in 1.0 1.5 2.0 2.5; do
    demb=$(mktemp --suffix=.dem)
    if python3 "$REPO/tools/make_synthetic_dem.py" --size "$W" --beta "$b" --seed 1 -o "$demb" >/dev/null 2>&1; then
      row=""
      for cap in "" 64 32 16 8 4 2; do
        out=$(DH_FLAT_PARTITION_REPLAY=1 "$BIN" "$demb" 0 "$SPEC8" $cap 2>/dev/null || true)
        t=$(grep -oE 'MPI-TREE-(MATCH|DIFFER)' <<<"$out" | head -1)
        v=$(grep -oE 'MPI-VOL-(MATCH|DIFFER)' <<<"$out" | head -1)
        if   [[ "$t" == *MATCH ]]; then c=M
        elif [[ "$v" == *MATCH ]]; then c=D
        else c='!'; fi
        row+="   $c"
      done
      printf "%-16s %s\n" "$b" "cap: $row"
    else printf "%-16s (generation failed)\n" "$b"; fi
    rm -f "$demb"
  done
else
  echo "   (skipped -- make_synthetic_dem.py unavailable; needs WTM_DIR)"
fi
