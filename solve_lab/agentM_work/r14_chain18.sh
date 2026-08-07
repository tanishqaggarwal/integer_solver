#!/bin/bash
# Start 2^18 ONLY after the 2^16 p80 run has finished (instruction: 2^16 first).
export PYTHONDONTWRITEBYTECODE=1
cd /home/user/integer_solver/solve_lab/agentM_work
until grep -q "^BEST" r14_enum16_p80.log; do sleep 30; done
echo "2^16 p80 complete; starting 2^18 at p80" > r14_chain18.log
exec python3 -u enumsub2.py 18 1e9 80 180 >> r14_enum18_p80.log 2>&1
