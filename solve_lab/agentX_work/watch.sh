#!/bin/sh
cd "$(dirname "$0")"
P=$(pgrep -x xmitm | head -1)
S0=$(awk '{print $14+$15}' /proc/$P/stat 2>/dev/null)
T0=$(date +%s)
sleep 300
S1=$(awk '{print $14+$15}' /proc/$P/stat 2>/dev/null)
T1=$(date +%s)
echo "pid=$P cpu_ticks_delta=$((S1-S0)) over $((T1-T0))s  => cores=$(echo "scale=2; ($S1-$S0)/100/($T1-$T0)"|bc)"
echo "i0 completions so far: $(grep -c i0= scanB.log)"
tail -3 scanB.log
