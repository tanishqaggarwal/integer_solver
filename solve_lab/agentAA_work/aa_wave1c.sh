#!/bin/bash
cd "$(dirname "$0")" || exit 1
./aa_run.sh tags_rem2.txt d 3 tbl/t4 tbl/bm4.bin
echo "=== PASS B (m<=7, a<=4 table) COMPLETE ==="
