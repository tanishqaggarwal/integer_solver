#!/bin/bash
# AH queue driver.  $1 = worker letter, $2 = budget seconds, rest = "n:seed" pairs.
export PYTHONDONTWRITEBYTECODE=1
cd /home/user/integer_solver/solve_lab/agentAH_work
WK=$1; BUD=$2; shift 2
for job in "$@"; do
  n=${job%%:*}; sd=${job##*:}
  tag="n${n}_s${sd}"
  if [ -f "meta_${tag}.json" ]; then echo "[$WK] skip $tag (done)"; continue; fi
  echo "[$WK] START $tag $(date +%H:%M:%S)"
  python3 ah_run.py "$tag" "$n" "$sd" 16 "$BUD"
  rc=$?
  echo "[$WK] EXIT $tag rc=$rc $(date +%H:%M:%S)"
done
echo "[$WK] QUEUE DRAINED $(date +%H:%M:%S)"
