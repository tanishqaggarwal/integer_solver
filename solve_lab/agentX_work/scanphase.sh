#!/bin/sh
cd "$(dirname "$0")"
DATA=$1; J=$2; REP=$3
B="0 5 10 17 25 38 128"
rm -f cpid_*
i=0; prev=""
for x in $B; do
  if [ -n "$prev" ]; then
    setsid nohup ./xrot scan $DATA $J xrot_tbl.bin xrot_bm.bin $REP $prev $x >/dev/null 2>&1 &
    echo $! > cpid_$i; i=$((i+1))
  fi
  prev=$x
done
sleep 1
for f in cpid_0 cpid_1 cpid_2 cpid_3 cpid_4 cpid_5; do p=$(cat $f); while kill -0 $p 2>/dev/null; do sleep 2; done; done
while pgrep -x xrot >/dev/null; do sleep 2; done
echo "SCAN PHASE DONE j=$J data=$DATA"
