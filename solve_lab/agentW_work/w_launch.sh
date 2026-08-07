#!/bin/bash
# launch $2.. detached, write the REAL python pid (not the wrapper's) to $1
cd /home/user/integer_solver/solve_lab/agentW_work
export PYTHONDONTWRITEBYTECODE=1
tag="$1"; shift
setsid nohup python3 "$@" > "${tag}.out" 2>&1 < /dev/null &
echo $! > "${tag}.pid"
