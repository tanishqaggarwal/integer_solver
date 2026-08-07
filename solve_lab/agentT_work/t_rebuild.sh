#!/bin/bash
# AUDIT T31 pre-step: the container restart wiped every *.pkl (they are gitignored, so the fresh
# checkout has none).  Rebuild the cache chain agentF_work + agentL_work need, INSIDE a private
# mirror under agentT_work/ so that nothing outside my directory is written.
set -e
export PYTHONDONTWRITEBYTECODE=1
MF=/home/user/integer_solver/solve_lab/agentT_work/mirror/F
ML=/home/user/integer_solver/solve_lab/agentT_work/mirror/L
cd $MF && echo "=== circ4"  && python3 -u circ4.py    | tail -3
cd $MF && echo "=== sched"  && python3 -u sched.py    | tail -3
cd $ML && echo "=== global(ors)"   && python3 -u global.py   | tail -3
cd $ML && echo "=== handles2"      && python3 -u handles2.py | tail -3
cd $ML && echo "=== buildall"      && python3 -u buildall.py | tail -4
cd $ML && echo "=== calib2"        && python3 -u calib2.py   | tail -4
cd $ML && echo "=== slopes"        && python3 -u slopes.py   | tail -3
echo "=== REBUILD DONE"
ls -la $MF/*.pkl $ML/*.pkl
