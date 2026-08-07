#!/bin/bash
cd "$(dirname "$0")" || exit 1
# Pass A: complete every remaining offset at m<=6 with the 89 MB a<=3 table (RAM-resident,
# immune to the fleet's page-cache pressure).
sed -e 's#runs/r_#runs6/r6_#' aa_run.sh > aa_run6.sh; chmod +x aa_run6.sh
./aa_run6.sh tags_rem.txt d 3 tbl/t3 tbl/bm3.bin
echo "=== PASS A (m<=6, a<=3 table) COMPLETE ==="
# Pass B: upgrade those same offsets to m<=7 with the big table, one at a time.
./aa_run.sh tags_rem.txt d 3 tbl/t4 tbl/bm4.bin
echo "=== PASS B (m<=7, a<=4 table) COMPLETE ==="
