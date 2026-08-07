#!/bin/bash
export PYTHONDONTWRITEBYTECODE=1
cd /home/user/integer_solver/solve_lab/agentU_work
rm -f u_pairs_*.pkl
for i in 0 1 2 3; do
  python3 -u u22_pairs.py $i 4 > u22_s$i.log 2>&1 &
  echo "shard $i pid $!"
done
wait
echo "ALL_SHARDS_DONE"
