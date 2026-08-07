#!/bin/bash
# wait (by PID, not by name) for the const tier to exit, then run the orbit tier
CP=$1
while kill -0 "$CP" 2>/dev/null; do sleep 20; done
echo "const pid $CP gone at $(date +%H:%M:%S)" >> chain.log
cd /home/user/integer_solver/solve_lab/agentAE_work
PYTHONDONTWRITEBYTECODE=1 python3 -u ae_families.py orbit R_orbit=48 log2max=27 > orbit48.log 2>&1
echo "orbit tier exit=$? at $(date +%H:%M:%S)" >> chain.log
