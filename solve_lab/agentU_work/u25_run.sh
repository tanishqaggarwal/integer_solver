#!/bin/bash
export PYTHONDONTWRITEBYTECODE=1
cd /home/user/integer_solver/solve_lab/agentU_work
rm -f u_exact_*.pkl
for i in 0 1 2 3; do python3 -u u25_exact.py $i 4 > u25_s$i.log 2>&1 & done
wait
echo ALL_DONE
