#!/bin/sh
cd "$(dirname "$0")"
export OMP_NUM_THREADS=8
export OMP_STACKSIZE=32M
for r in "0 8" "8 16" "16 32" "32 64" "64 128" "128 256"; do
  set -- $r
  ./xmitm scan data_real.txt 5 tbl4s.bin bm4.bin rep_real.txt $1 $2
  echo "=== chunk [$1,$2) COMPLETE ==="
done
echo "=== PASS B COMPLETE (|S| <= 9 exhausted) ==="
