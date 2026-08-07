#!/bin/bash
# Agent Y -- the real sweep against the COMPLEMENT target T'.
# Table = agent X's |A| in 1..4 key table (target-independent, revalidated by ytable_check.py).
# Coverage: scan size b + table size a covers |S'| = b+a for a in 1..4.
#   size 2 -> |S'| 3..6 ;  size 3 -> 4..7 ;  size 4 -> 5..8 ;  size 5 -> 6..9
#   |S'| <= 4 is closed by yedge.py (direct table probe of T').
cd "$(dirname "$0")"
TBL=../agentX_work/tbl4s.bin
BM=../agentX_work/bm4.bin
export OMP_NUM_THREADS=4
for SZ in 2 3 4; do
  ./ymitm scan data_comp.txt $SZ $TBL $BM rep_comp.txt >>yrun_$SZ.log 2>&1
  echo "finished size $SZ at $(date -u +%H:%M:%S)" >> yrun.status
done
./ymitm scan data_comp.txt 5 $TBL $BM rep_comp.txt >>yrun_5.log 2>&1
echo "finished size 5 at $(date -u +%H:%M:%S)" >> yrun.status
echo "ALLDONE" >> yrun.status
