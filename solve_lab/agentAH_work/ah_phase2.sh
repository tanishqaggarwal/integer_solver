#!/bin/bash
# Phase 2: bounded-representation variant (ah_joint), 3 GB cap, exact fast root enumeration.
# Validation points first (must reproduce the verbatim runs byte-for-byte), then the grid.
cd /home/user/integer_solver/solve_lab/agentAH_work
export PYTHONDONTWRITEBYTECODE=1 AH_FASTROOTS=1 AH_BIGROOT=1 AH_MEMCAP_GB=3 AH_TAGSUF=g AH_OUTER=60
nohup ./ah_driver.sh A 3600 32:101 8:101 96:101 128:101 160:101 192:101 24:101 48:101 64:101 \
      96:202 128:202 160:202 192:202 24:202 48:202 64:202 32:202 16:202 8:202 4:202 2:202 1:202 \
      > drvAg.log 2>&1 &
echo "PID_Ag=$!"
sleep 3
nohup ./ah_driver.sh B 3600 255:101 252:101 248:101 240:101 224:101 255:202 252:202 248:202 240:202 224:202 \
      255:303 252:303 248:303 240:303 224:303 192:303 160:303 128:303 96:303 \
      > drvBg.log 2>&1 &
echo "PID_Bg=$!"
