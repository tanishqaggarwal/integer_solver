#!/bin/bash
cd "$(dirname "$0")" || exit 1
while read -r nm s0; do
  for s in 0 1 2 3 4 5 6 7; do
    grep -q "SHARDOK$s" "runs/q_${nm}.done" 2>/dev/null && continue
    OMP_NUM_THREADS=4 ./aa_signed scan "data/q_$nm.txt" 4 "tbl/v$s/t4" tbl/bm4.bin \
        "runs/q_$nm.txt" "$s0" $((s0+1)) 2> "shardlogs/q_$nm.s$s.err"
    rc=$?
    [ $rc -eq 0 ] && echo "SHARDOK$s" >> "runs/q_${nm}.done" || echo "FAILED $nm shard=$s exit=$rc"
  done
  echo "$nm engineDONE=$(grep -c DONE runs/q_$nm.txt) hits=$(grep -c HIT runs/q_$nm.txt)"
done < plant8.list
echo "PLANT8_END"
