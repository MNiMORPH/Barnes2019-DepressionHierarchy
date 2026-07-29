#!/usr/bin/env bash
# check_split_invariance.sh -- the split-invariance acceptance check (see dephier_stitch STITCH-SIG).
# Builds the tiled DepressionHierarchy at two-or-more different tilings and asserts the tree is IDENTICAL
# across them -- "identical DH regardless of how we split", the real acceptance metric (does NOT reference
# serial). Exits 0 and prints SPLIT-INVARIANT iff every tiling yields the same STITCH-SIG; else exits 1 and
# prints SPLIT-VARIANT with the distinct signatures.
#
# Usage: check_split_invariance.sh <dephier_stitch.exe> <dem> <ocean_level> <splitsA> <splitsB> [<splitsC> ...]
#   e.g. check_split_invariance.sh ./dephier_stitch.exe kerry_test7.dem 0 2 3 5 7
set -u
exe=$1; dem=$2; ocean=$3; shift 3
[ $# -ge 2 ] || { echo "need >=2 split configs to compare" >&2; exit 2; }

sigs=""
for splits in "$@"; do
  s=$("$exe" "$dem" "$ocean" "$splits" 2>/dev/null | grep -aoE "STITCH-SIG .* sig=[0-9]+" | grep -aoE "sig=[0-9]+")
  [ -n "$s" ] || { echo "no STITCH-SIG for splits=$splits (build failed/aborted?)" >&2; exit 2; }
  sigs="$sigs $splits:$s"
done

uniq=$(for kv in $sigs; do echo "${kv#*:}"; done | sort -u | wc -l)
if [ "$uniq" -eq 1 ]; then
  echo "SPLIT-INVARIANT $dem  [$sigs ]"
  exit 0
else
  echo "SPLIT-VARIANT $dem  ($uniq distinct trees) [$sigs ]"
  exit 1
fi
