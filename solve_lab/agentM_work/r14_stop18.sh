#!/bin/bash
# Coordinator: let 2^18 run to |W| = 8 and then stop it (39,026 never occurs above |W| = 7).
cd /home/user/integer_solver/solve_lab/agentM_work
until grep -q "^  |W|= 8         COMPLETE" r14_enum18_p80.log; do
  kill -0 28848 2>/dev/null || { echo "child already gone" >> r14_stop18.log; exit 0; }
  sleep 30
done
echo "|W|=8 complete -> stopping 2^18 (pid 28848)" >> r14_stop18.log
kill 28848 2>/dev/null; sleep 3; kill 28845 2>/dev/null
echo "stopped" >> r14_stop18.log
