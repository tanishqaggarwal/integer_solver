#!/bin/sh
cd "$(dirname "$0")"
./xmitm scan data_real.txt 2 tbl4s.bin bm4.bin rep_real.txt
./xmitm scan data_real.txt 3 tbl4s.bin bm4.bin rep_real.txt
./xmitm scan data_real.txt 4 tbl4s.bin bm4.bin rep_real.txt
echo "=== PASS A COMPLETE (|S| <= 8 exhausted) ==="
./xmitm scan data_real.txt 5 tbl4s.bin bm4.bin rep_real.txt
echo "=== PASS B COMPLETE (|S| <= 9 exhausted) ==="
