#!/bin/bash
# Wall-clock and peak-RSS benchmark: R pipeline vs the C executable.
set -u
BIN=../physmerge
DIR=${DIR:-/private/tmp/physmerge_bench}
REPS=${REPS:-3}
mkdir -p "$DIR"
run() {   # run() label cmd... -> "label  median_s  peak_MB"
  local label="$1"; shift
  local best="" ; local times=()
  for r in $(seq 1 "$REPS"); do
    out=$( { /usr/bin/time -l "$@" >/dev/null; } 2>&1 )
    t=$(echo "$out" | awk '/real/{print $1; exit}')
    m=$(echo "$out" | awk '/maximum resident set size/{print $1; exit}')
    times+=("$t"); best="$m"
  done
  med=$(printf '%s\n' "${times[@]}" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  printf '%-28s %8.2f s %10.1f MB\n' "$label" "$med" "$(echo "$best/1048576" | bc -l)"
}
for N in "$@"; do
  F="$DIR/gwas_${N}.tsv"
  [ -f "$F" ] || Rscript gen_bench.R "$N" "$F"
  SZ=$(du -m "$F" | cut -f1)
  echo "=== n=$N rows (${SZ} MB) ==="
  run "R pipeline"  Rscript run_r.R "$F" "$DIR/r_ids_${N}.txt"
  run "C executable" $BIN --input "$F" --format plink2 --sig-th 5e-8 --window 500000 \
        --reset-on any --quiet --out "$DIR/c_blocks_${N}.tsv" --snp-list "$DIR/c_ids_${N}.txt"
  if diff -q <(sort "$DIR/r_ids_${N}.txt") <(sort "$DIR/c_ids_${N}.txt") >/dev/null; then
    echo "  lead-SNP lists identical ($(wc -l < "$DIR/c_ids_${N}.txt" | tr -d ' ') blocks)"
  else
    echo "  !! lead-SNP lists DIFFER"
  fi
  echo
done
