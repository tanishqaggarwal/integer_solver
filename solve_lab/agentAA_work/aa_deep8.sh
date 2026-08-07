#!/bin/bash
# Agent AA -- c0 at m <= 8:  a<=4 table x b=4 scan, the only split that reaches m=8.
#
# Structure: OUTER loop = s0 chunk (16 balanced ranges, chunks4.json), INNER loop = table shard
# (8 single-shard passes whose union is the whole a<=4 table).  All eight shards of a chunk are
# finished before the next chunk starts, so "chunks 0..j complete" is an EXACT fraction of the
# m=8 candidate space, not a vague partial.
#
# Evidence discipline: engine stderr to a per-unit log; nothing but the engine writes the
# evidence file; exit code tested; resume keys on the engine's own DONE line carrying the exact
# expected candidate count for that [lo,hi) range.  A crashed unit looks incomplete, which is
# the safe direction.
cd "$(dirname "$0")" || exit 1
mkdir -p runs8 shardlogs8
NCH=$(python3 -c "import json;print(len(json.load(open('chunks4.json'))))")
for ((j=0; j<NCH; j++)); do
  read -r LO HI NEXP < <(python3 -c "
import json
from math import comb
ch=json.load(open('chunks4.json'))[$j]
n=sum(comb(255-(s>>1),3)*8 for s in range(ch[0],ch[1]))
print(ch[0],ch[1],n)")
  for s in 0 1 2 3 4 5 6 7; do
    ev="runs8/c0m8.c${j}.s${s}.txt"
    grep -q "range=\[$LO,$HI) n=$NEXP " "$ev" 2>/dev/null && continue
    OMP_NUM_THREADS=4 ./aa_signed scan data/d_c0.txt 4 "tbl/v$s/t4" tbl/bm4.bin \
        "$ev" "$LO" "$HI" 2> "shardlogs8/c${j}.s${s}.err"
    rc=$?
    if [ $rc -ne 0 ]; then echo "FAILED chunk=$j shard=$s exit=$rc"; continue; fi
    grep -q "range=\[$LO,$HI) n=$NEXP " "$ev" || echo "NOCOUNT chunk=$j shard=$s"
  done
  ok=0
  for s in 0 1 2 3 4 5 6 7; do
    grep -q "range=\[$LO,$HI) n=$NEXP " "runs8/c0m8.c${j}.s${s}.txt" 2>/dev/null && ok=$((ok+1))
  done
  echo "chunk=$j range=[$LO,$HI) n_each=$NEXP shards_with_exact_count=$ok/8 hits=$(cat runs8/c0m8.c${j}.s*.txt 2>/dev/null | grep -c HIT)"
done
echo "DEEP8_END"
