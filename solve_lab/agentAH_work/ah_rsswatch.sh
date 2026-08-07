#!/bin/bash
# Watchdog for the two ah_run processes that started BEFORE the RLIMIT_AS guard existed.
# Kills a process whose RSS exceeds the cap, so the kernel OOM killer never has to choose
# between my job and another agent's.  Identifies processes by PID only.
CAP_KB=3200000
for P in "$@"; do :; done
while true; do
  live=0
  for P in "$@"; do
    kill -0 "$P" 2>/dev/null || continue
    live=1
    rss=$(awk '/VmRSS/{print $2}' /proc/$P/status 2>/dev/null)
    [ -z "$rss" ] && continue
    if [ "$rss" -gt "$CAP_KB" ]; then
      echo "$(date +%H:%M:%S) RSS-WATCHDOG killing pid=$P rss_kB=$rss > $CAP_KB"
      kill -9 "$P"
    fi
  done
  [ "$live" -eq 0 ] && break
  sleep 5
done
echo "$(date +%H:%M:%S) watchdog: all watched pids gone"
