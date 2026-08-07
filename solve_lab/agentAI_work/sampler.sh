#!/bin/sh
echo $$ > /home/user/integer_solver/solve_lab/agentAI_work/sampler.pid
cd /home/user/integer_solver/solve_lab/agentAI_work
X=/home/user/integer_solver/solve_lab/agentX_work
while :; do
  TS=$(date -u +%H:%M:%S)
  AVAIL=$(df -m /home/user | tail -1 | awk '{print $4}')
  AA=$(du -sm /home/user/integer_solver/solve_lab/agentAA_work 2>/dev/null | cut -f1)
  XW=$(du -sm $X 2>/dev/null | cut -f1)
  ROT=$(grep -c ' DONE$' $X/rotdone.txt 2>/dev/null)
  SH_=$(grep -c '^DONE rot=' $X/rrep_real.txt 2>/dev/null)
  HIT=$(grep -c '^HIT' $X/rrep_real.txt 2>/dev/null)
  kill -0 30892 2>/dev/null && A=alive || A=DEAD
  if kill -0 6881 2>/dev/null; then B=alive; BR=$(awk '/VmRSS/{print $2/1024}' /proc/6881/status 2>/dev/null | cut -d. -f1); else B=DEAD; BR=0; fi
  MAV=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
  LOAD=$(cut -d' ' -f1 /proc/loadavg)
  echo "$TS availMB=$AVAIL aaMB=$AA xMB=$XW rot=$ROT shards=$SH_ hits=$HIT A(30892)=$A B(6881)=$B dregRSSmb=$BR memAvailMB=$MAV load=$LOAD" >> sampler.log
  # alarms -> separate file so they are easy to spot
  [ "$AVAIL" -lt 6144 ] && echo "$TS ALARM disk availMB=$AVAIL < 6144" >> ALARMS.log
  [ "$MAV" -lt 1024 ] && echo "$TS ALARM memAvailable=${MAV}MB < 1024" >> ALARMS.log
  [ "$A" = DEAD ] && echo "$TS ALARM rotall(30892) DEAD rot=$ROT" >> ALARMS.log
  [ "$B" = DEAD ] && echo "$TS ALARM dreg(6881) DEAD" >> ALARMS.log
  [ "$HIT" -gt 0 ] && echo "$TS ALARM HIT lines present in rrep_real.txt: $HIT" >> ALARMS.log
  sleep 120
done
