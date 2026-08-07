#!/bin/sh
# one rotation of the 128-rotation splitting sweep.  usage: rotone.sh <j> <repfile>
cd "$(dirname "$0")"
J=$1; REP=$2
B="0 5 10 17 25 38 128"
waitpids(){ for p in $(cat "$@" 2>/dev/null); do while kill -0 "$p" 2>/dev/null; do sleep 2; done; done; }
# --- build phase: 6 sessions ---
rm -f rt_*.bin rtpid_*
i=0; prev=""
for x in $B; do
  if [ -n "$prev" ]; then
    setsid nohup ./xrot build data_real.txt $J $prev $x rt_$i.bin >/dev/null 2>&1 &
    echo $! > rtpid_$i; i=$((i+1))
  fi
  prev=$x
done
sleep 1
waitpids rtpid_0 rtpid_1 rtpid_2 rtpid_3 rtpid_4 rtpid_5
while pgrep -x xrot >/dev/null; do sleep 2; done
# --- sort phase: 6 sessions ---
rm -f spid_*
for i in 0 1 2 3 4 5; do
  setsid env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$PWD/pylib nohup python3 sortchunk.py rt_$i.bin >/dev/null 2>&1 &
  echo $! > spid_$i
done
sleep 1
waitpids spid_0 spid_1 spid_2 spid_3 spid_4 spid_5
# --- merge + bitmap ---
./xrot merge rt_all.bin rt_0.bin rt_1.bin rt_2.bin rt_3.bin rt_4.bin rt_5.bin 2>/dev/null
rm -f rt_0.bin rt_1.bin rt_2.bin rt_3.bin rt_4.bin rt_5.bin
./xrot bitmap rt_all.bin rt_bm.bin 2>/dev/null
# --- scan phase: 6 sessions ---
rm -f cpid_*
i=0; prev=""
for x in $B; do
  if [ -n "$prev" ]; then
    setsid nohup ./xrot scan data_real.txt $J rt_all.bin rt_bm.bin $REP $prev $x >/dev/null 2>&1 &
    echo $! > cpid_$i; i=$((i+1))
  fi
  prev=$x
done
sleep 1
waitpids cpid_0 cpid_1 cpid_2 cpid_3 cpid_4 cpid_5
while pgrep -x xrot >/dev/null; do sleep 2; done
rm -f rt_all.bin rt_bm.bin
echo "ROTATION $J COMPLETE"
