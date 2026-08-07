#!/bin/bash
# m<=7 via 8 single-shard passes: working set 1.4 GB shard + 0.5 GB bitmap, always resident.
# Union of the 8 passes == the full a<=4 table, since shard = key>>61 partitions the keyspace.
cd "$(dirname "$0")" || exit 1
while read -r tag; do
  [ -z "$tag" ] && continue
  out=runs/rs_${tag}.txt
  for b in 1 2 3; do
    for s in 0 1 2 3 4 5 6 7; do
      grep -q "sz=$b .*SHARD$s\$" "$out" 2>/dev/null && continue
      OMP_NUM_THREADS=4 ./aa_signed scan "data/d_${tag}.txt" "$b" "tbl/v$s/t4" tbl/bm4.bin "$out" >/dev/null 2>&1
      echo "SHARD$s" >> "$out"
    done
  done
  echo "SHARDED_DONE $tag hits=$(grep -c HIT "$out" 2>/dev/null) lines=$(grep -c DONE "$out")"
done < "$1"
echo "SHARDED_COMPLETE"
