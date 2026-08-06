#!/bin/bash
# chunked driver: 4 workers x 420 candidates per round, sequential rounds
cd /home/user/integer_solver/solve_lab
for start in 1680 3360 5040; do
  for k in 0 1 2 3; do
    s=$((start + k*420)); e=$((s+420))
    timeout 900 python3 s10/bn_key4.py $s $e > s10/bn_key4_$s.log 2>&1 &
  done
  wait
  echo "round $start complete"
done
echo ALLDONE
