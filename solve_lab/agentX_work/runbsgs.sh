#!/bin/sh
cd "$(dirname "$0")"
cat babys.bin > /dev/null
echo "=== planted test (expect CAND i=5, k0=336778887) ==="
nice -n 19 ./xbsgs giant plant_bsgs.txt babys.bin 1 2>/dev/null
echo "=== real target, +T ==="
nice -n 19 ./xbsgs giant data_real.txt babys.bin 1 2>/dev/null
echo "=== real target, -T ==="
nice -n 19 ./xbsgs giant data_real.txt babys.bin -1 2>/dev/null
echo "=== BSGS ALL DONE ==="
