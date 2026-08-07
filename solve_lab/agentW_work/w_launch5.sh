#!/bin/bash
cd /home/user/integer_solver/solve_lab/agentW_work
export PYTHONDONTWRITEBYTECODE=1 WFILES=w_cocirc3_raw_s1_3.json WOUT=w_close3_s13.json
setsid nohup python3 w_close3.py > w_close3_s13.out 2>&1 < /dev/null &
echo $! > w_close3_s13.pid
