#!/bin/bash
# m<=7 via 8 single-shard passes.  Working set = 1.4 GB shard + 0.5 GB bitmap, always resident,
# so this is immune to the page-cache pressure that made the monolithic 11.3 GB table disk-bound.
# The union of the 8 passes IS the full a<=4 table, because shard = key>>61 partitions the
# 64-bit key space and every pass enumerates the identical candidate set.
#
# EVIDENCE DISCIPLINE (coordinator's rule, after agent Y's yorbit.status and agent AI's audit):
#   * the engine's stderr goes to a per-shard log, never to /dev/null;
#   * NOTHING is written to the evidence file by this shell -- only the engine writes there;
#   * the exit code is tested, and a non-zero exit is recorded as FAILED and re-tried on resume;
#   * the resume guard keys on the ENGINE's own "DONE ... n=<count>" line for that exact shard,
#     checked against the closed form C(256,b)*2^b -- never on a marker this script wrote.
# A shard that crashes therefore looks incomplete on resume, which is the direction that is safe.
cd "$(dirname "$0")" || exit 1
PRE=${2:-d}
exp_n () { python3 -c "from math import comb;print(comb(256,$1)*2**$1)"; }
declare -A EXP
for b in 1 2 3; do EXP[$b]=$(exp_n "$b"); done
mkdir -p runs shardlogs
fail=0
while read -r tag; do
  [ -z "$tag" ] && continue
  for b in 1 2 3; do
    for s in 0 1 2 3 4 5 6 7; do
      ev="runs/rs_${PRE}${tag}.s${s}.txt"
      # resume guard: engine evidence only
      if grep -q "sz=$b .* n=${EXP[$b]} " "$ev" 2>/dev/null; then continue; fi
      OMP_NUM_THREADS=4 ./aa_signed scan "data/${PRE}_${tag}.txt" "$b" "tbl/v$s/t4" tbl/bm4.bin \
          "$ev" 2> "shardlogs/${PRE}${tag}.b${b}.s${s}.err"
      rc=$?
      if [ $rc -ne 0 ]; then
        echo "FAILED tag=$tag b=$b shard=$s exit=$rc  (see shardlogs/${PRE}${tag}.b${b}.s${s}.err)"
        fail=$((fail+1)); continue
      fi
      # verify the engine really produced the exhaustive count before moving on
      if ! grep -q "sz=$b .* n=${EXP[$b]} " "$ev"; then
        echo "NOCOUNT tag=$tag b=$b shard=$s exit=0 but no exhaustive DONE line"
        fail=$((fail+1))
      fi
    done
  done
  # summary line derived from evidence, not from having reached this point in the loop
  ok=0
  for b in 1 2 3; do for s in 0 1 2 3 4 5 6 7; do
    grep -q "sz=$b .* n=${EXP[$b]} " "runs/rs_${PRE}${tag}.s${s}.txt" 2>/dev/null && ok=$((ok+1))
  done; done
  echo "tag=$tag exhaustive_engine_DONE_lines=$ok/24 hits=$(cat runs/rs_${PRE}${tag}.s*.txt 2>/dev/null | grep -c HIT)"
done < "$1"
echo "SHARDED_END failures=$fail"
