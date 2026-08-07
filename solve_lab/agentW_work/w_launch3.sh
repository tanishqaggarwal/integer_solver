#!/bin/bash
cd /home/user/integer_solver/solve_lab/agentW_work
export PYTHONDONTWRITEBYTECODE=1 WSMIN=4 WSMAX=6
setsid nohup python3 w_cocirc3.py > w_cocirc3.out 2>&1 < /dev/null &
echo $! > w_cocirc3.pid
