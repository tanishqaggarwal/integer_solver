#!/bin/bash
# Agent AA -- sequential offset sweep.  ONE compute process at a time.
#   $1 = tag list file (one offset tag per line)
#   $2 = prefix for data files: 'd' (real) or 'p' (plant)
#   $3 = max scan size b
#   $4 = table prefix
#   $5 = bitmap
cd "$(dirname "$0")" || exit 1
LIST=$1; PRE=$2; BMAX=$3; TBL=$4; BM=$5
mkdir -p runs
while read -r tag; do
  [ -z "$tag" ] && continue
  out=runs/r_${PRE}_${tag}.txt
  for b in $(seq 1 "$BMAX"); do
    grep -q "sz=$b " "$out" 2>/dev/null && continue
    OMP_NUM_THREADS=4 ./aa_signed scan "data/${PRE}_${tag}.txt" "$b" "$TBL" "$BM" "$out" 2>&1 \
      | tail -1
  done
  echo "OFFSET_DONE $tag hits=$(grep -c HIT "$out" 2>/dev/null)"
done < "$LIST"
echo "SWEEP_COMPLETE"
