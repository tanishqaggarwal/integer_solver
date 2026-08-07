#!/bin/bash
# Wait for the two orphaned verbatim runs to exit, then start phase 2 automatically.
cd /home/user/integer_solver/solve_lab/agentAH_work
while kill -0 9740 2>/dev/null || kill -0 13873 2>/dev/null; do sleep 20; done
echo "$(date +%H:%M:%S) both verbatim high-|S| runs finished; starting phase 2"
./ah_phase2.sh
