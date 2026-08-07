#!/bin/bash
cd "$(dirname "$0")" || exit 1
echo "### VALIDATION of the 8-shard union path (must find planted m=5 answers)"
./aa_shard.sh tags_plantshard.txt p
echo "### REMAINING OFFSETS at m<=7 via shard passes"
./aa_shard.sh tags_rem3.txt d
echo "### WAVE3 END"
