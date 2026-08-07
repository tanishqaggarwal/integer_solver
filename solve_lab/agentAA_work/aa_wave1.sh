#!/bin/bash
cd "$(dirname "$0")" || exit 1
echo "=== PLANT VALIDATION against a<=4 table ==="
./aa_run.sh tags_plants.txt p 3 tbl/t4 tbl/bm4.bin
echo "=== REAL OFFSET SWEEP m<=7 ==="
./aa_run.sh tags_all.txt d 3 tbl/t4 tbl/bm4.bin
echo "=== WAVE1 ALL DONE ==="
