#!/bin/sh
# Custodial sampler -- agent AI.  Job A (X's rotational sweep) only.
#
# DESIGN RULE, learned the hard way (see ALARMS_QUARANTINE_README.md):
#   Absence of a process is NOT evidence of failure.  `kill -0` returning false
#   says only "the process is gone" -- equally consistent with success and with
#   failure.  So on the PID going absent this script READS THE ARTEFACT and only
#   then names COMPLETED or DIED, quoting the terminal line either way.
#
# PID HANDLE (agent AA's rule): a recorded PID for a job that respawns its worker
#   is stale by construction; the durable handle is the PARENT.  30892 is
#   rotall.sh, the parent that respawns the six ./xrot workers every rotation.
#   The workers are deliberately NOT tracked by PID.
#
# Alarms are LATCHED -- each condition fires at most once, so one inference is
#   never restated until it looks like many.
#
# Job B (dreg, PID 6881) is NOT watched: it completed at 22:29:50 and is recorded.

D=/home/user/integer_solver/solve_lab/agentAI_work
X=/home/user/integer_solver/solve_lab/agentX_work
echo $$ > $D/sampler.pid
cd $D

APID=30892                 # rotall.sh -- the PARENT, the durable handle
TERM_A='ALL 128 ROTATIONS ATTEMPTED'
L_DISK=0; L_MEM=0; L_HIT=0

alarm(){ echo "$(date -u +%H:%M:%S) $*" >> $D/ALARMS.log; }

while :; do
  TS=$(date -u +%H:%M:%S)
  AVAIL=$(df -m /home/user | tail -1 | awk '{print $4}')
  XW=$(du -sm $X 2>/dev/null | cut -f1)
  ROT=$(grep -c ' DONE$' $X/rotdone.txt 2>/dev/null)
  INC=$(grep -c 'INCOMPLETE' $X/rotdone.txt 2>/dev/null)
  SH_=$(grep -c '^DONE rot=' $X/rrep_real.txt 2>/dev/null)
  HIT=$(grep -c '^HIT' $X/rrep_real.txt 2>/dev/null)
  MAV=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
  LOAD=$(cut -d' ' -f1 /proc/loadavg)

  if kill -0 $APID 2>/dev/null; then
    # confirm the recorded PID still refers to rotall.sh (guards against PID reuse;
    # this is confirming a KNOWN pid, not discovering one by matching)
    case "$(tr '\0' ' ' < /proc/$APID/cmdline 2>/dev/null)" in
      *rotall.sh*) A=alive ;;
      *)           A=REUSED ;;
    esac
  else
    A=absent
  fi

  echo "$TS availMB=$AVAIL xMB=$XW rot=$ROT inc=$INC shards=$SH_ hits=$HIT A($APID)=$A memAvailMB=$MAV load=$LOAD" >> $D/sampler.log

  # ---- Job A terminal handling: READ THE ARTEFACT before naming the outcome ----
  if [ "$A" != alive ]; then
    if grep -q "$TERM_A" $X/rotdone.txt 2>/dev/null; then
      alarm "JOB A COMPLETED (pid $APID absent AND terminal line present)"
      alarm "  terminal line: $(grep "$TERM_A" $X/rotdone.txt | tail -1)"
      alarm "  rotations marked DONE: $ROT/128   INCOMPLETE: $INC   HIT lines: $HIT"
      alarm "  ACTION: run verify_rot.py; |S|<=10 claimable ONLY if 128 DONE and it exits 0"
    else
      alarm "JOB A DIED (pid $APID absent AND terminal line '$TERM_A' ABSENT)"
      alarm "  last rotdone.txt line: $(tail -1 $X/rotdone.txt 2>/dev/null)"
      alarm "  progress at death: rot=$ROT/128 INCOMPLETE=$INC shards=$SH_ hits=$HIT"
      alarm "  ACTION: sweep is restartable and skips rotations already DONE -- restart it,"
      alarm "          and record in CUSTODY.md that it died, when, and why, with evidence."
    fi
    exit 0   # stop watching a job that has reached a terminal state
  fi

  # ---- latched resource / result alarms ----
  if [ "$AVAIL" -lt 6144 ] && [ $L_DISK -eq 0 ]; then
    alarm "disk availMB=$AVAIL < 6144 (trough of X's per-rotation sawtooth; recommend only, never delete)"; L_DISK=1
  fi
  [ "$AVAIL" -ge 7168 ] && L_DISK=0          # re-arm once clearly recovered
  if [ "$MAV" -lt 1024 ] && [ $L_MEM -eq 0 ]; then
    alarm "memAvailable=${MAV}MB < 1024"; L_MEM=1
  fi
  [ "$MAV" -ge 2048 ] && L_MEM=0
  if [ "$HIT" -gt 0 ] && [ $L_HIT -eq 0 ]; then
    alarm "HIT lines present in rrep_real.txt: $HIT -- a candidate S may have been found; verify before claiming"; L_HIT=1
  fi

  sleep 120
done
