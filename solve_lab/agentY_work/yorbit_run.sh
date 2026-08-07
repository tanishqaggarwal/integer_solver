#!/bin/bash
# Agent Y -- the endomorphism orbit sweep.
#   ./yorbit_run.sh 4      -> scan sizes 2,3,4  (covers |S| <= 8)   ~25 s per target
#   ./yorbit_run.sh 5      -> scan size 5 only  (extends to |S| <= 9) ~18 min per target
# T itself was swept to |S| <= 9 by agent X; c_T is the complement target already swept here.
cd "$(dirname "$0")"
TBL=../agentX_work/tbl4s.bin
BM=../agentX_work/bm4.bin
export OMP_NUM_THREADS=4
MAX=${1:-4}
for NM in negT lamT neglamT lam2T neglam2T c_negT c_lamT c_neglamT c_lam2T c_neglam2T; do
  if [ "$MAX" = "5" ]; then SIZES="5"; else SIZES="2 3 4"; fi
  for SZ in $SIZES; do
    ./ymitm scan data_$NM.txt $SZ $TBL $BM rep_orbit_$NM.txt >>yorbit_$NM.log 2>&1
  done
  echo "$NM done sizes=$SIZES at $(date -u +%H:%M:%S)" >> yorbit.status
done
echo "ORBITDONE_$MAX" >> yorbit.status
