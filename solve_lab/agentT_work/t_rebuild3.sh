#!/bin/bash
# T35 (post second restart): the FULL chain.  t_rebuild.sh omitted ortree2.py and handles.py,
# which handles2/buildall depend on -- ran by hand last time.  This is the complete script.
set -e
export PYTHONDONTWRITEBYTECODE=1
MF=/home/user/integer_solver/solve_lab/agentT_work/mirror/F
ML=/home/user/integer_solver/solve_lab/agentT_work/mirror/L
cd $MF; [ -f circ4.pkl ] || { echo "=== circ4"; python3 -u circ4.py | tail -3; }
cd $MF; [ -f sched.pkl ] || { echo "=== sched"; python3 -u sched.py | tail -3; }
cd $ML; [ -f ors.pkl ]   || { echo "=== global"; python3 -u global.py | tail -3; }
cd $ML; echo "=== ortree2" ; python3 -u ortree2.py  2>&1 | tail -4
cd $ML; echo "=== handles" ; python3 -u handles.py  2>&1 | tail -4
cd $ML; echo "=== handles2"; python3 -u handles2.py 2>&1 | tail -4
cd $ML; echo "=== buildall"; python3 -u buildall.py 2>&1 | tail -4
cd $ML; echo "=== calib2"  ; python3 -u calib2.py   2>&1 | tail -5
cd $ML; echo "=== slopes"  ; python3 -u slopes.py   2>&1 | tail -3
echo "=== REBUILD DONE"; ls -la $MF/*.pkl $ML/*.pkl
