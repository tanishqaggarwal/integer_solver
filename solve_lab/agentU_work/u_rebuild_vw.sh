#!/bin/bash
set -e
export PYTHONDONTWRITEBYTECODE=1
cd /home/user/integer_solver/solve_lab/agentU_work
for s in v1_parse.py v3_defs.py v5_chain.py v8b_supp.py v9b_theorem.py w1_zfactors.py w2_wire.py w3_crt.py; do
  echo "=== $s"; python3 -u $s 2>&1 | tail -6
done
echo "=== VW_REBUILD_DONE"
