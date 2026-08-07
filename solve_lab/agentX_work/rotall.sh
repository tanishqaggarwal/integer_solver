#!/bin/sh
cd "$(dirname "$0")"
for J in $(seq 0 127); do
  grep -q "^ROT $J DONE" rotdone.txt 2>/dev/null && continue
  ./rotone.sh $J rrep_real.txt >/dev/null 2>&1
  n=$(grep -c "rot=$J " rrep_real.txt)
  if [ "$n" = "6" ]; then echo "ROT $J DONE" >> rotdone.txt; else echo "ROT $J INCOMPLETE ($n/6)" >> rotdone.txt; fi
done
echo "ALL 128 ROTATIONS ATTEMPTED" >> rotdone.txt
