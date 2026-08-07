#!/bin/bash
cd /home/user/integer_solver
for f in solve_lab/agentM_work/xc14_*.json; do
  exp=$(basename $f | cut -d_ -f2)
  line=$(timeout 900 python3 solve_lab/checker.py "$f" 2>&1 | grep "satisfied")
  got=$(echo "$line" | sed -E 's/.*satisfied ([0-9]+)\/.*/\1/')
  if [ "$exp" == "$got" ]; then v=AGREE; else v="DISAGREE"; fi
  echo "$v  engine=$exp  checker=$got  $(basename $f)"
done
