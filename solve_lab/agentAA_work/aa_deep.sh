#!/bin/bash
# b=4 deep sweep (m<=8) on a priority list, chunked by s0 so partial coverage is an exact
# fraction.  ONE process at a time.  $1 = file of tags.
cd "$(dirname "$0")" || exit 1
CH=$(python3 -c "import json;print(' '.join('%d:%d'%tuple(x) for x in json.load(open('chunks4.json'))))")
while read -r tag; do
  [ -z "$tag" ] && continue
  out=runs/r4_${tag}.txt
  for c in $CH; do
    lo=${c%%:*}; hi=${c##*:}
    grep -q "range=\[$lo,$hi)" "$out" 2>/dev/null && continue
    OMP_NUM_THREADS=4 ./aa_signed scan "data/d_${tag}.txt" 4 tbl/t4 tbl/bm4.bin "$out" "$lo" "$hi" 2>&1 | tail -1
  done
  echo "DEEP_DONE $tag hits=$(grep -c HIT "$out" 2>/dev/null)"
done < "$1"
echo "DEEP_COMPLETE"
