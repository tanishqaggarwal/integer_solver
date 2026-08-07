#!/bin/bash
# Records the R=60 verdict when the engine exits.  Waits by PID, never by name.
# Writes R60_VERDICT.txt; if a candidate appears it is verified by the INDEPENDENT
# implementation (ae_verify.py) before anything is written.
cd /home/user/integer_solver/solve_lab/agentAE_work
EPID=$1
while kill -0 "$EPID" 2>/dev/null; do sleep 30; done
{
  echo "engine pid $EPID exited at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if grep -q '^CAND' head60b.out 2>/dev/null; then
      echo "STATUS: CANDIDATE PRESENT -- NOT A RESULT UNTIL VERIFIED"
      grep '^CAND' head60b.out
      echo "run:  python3 ae_lib.py-side reconstruction, then python3 ae_verify.py <k>"
  elif grep -q '^DONE' head60b.out 2>/dev/null; then
      D=$(grep '^DONE' head60b.out)
      echo "$D"
      J=$(sed 's/.*jumps=\([0-9]*\).*/\1/' <<<"$D")
      DP=$(sed 's/.*dps=\([0-9]*\).*/\1/' <<<"$D")
      DX=$(sed 's/.*dpexp=\([0-9.]*\).*/\1/' <<<"$D")
      C=$(sed 's/.*cands=\([0-9]*\).*/\1/' <<<"$D")
      python3 - "$J" "$DP" "$DX" "$C" <<'PY'
import sys, math
J=int(sys.argv[1]); DP=int(sys.argv[2]); DX=float(sys.argv[3]); C=int(sys.argv[4])
print("  jumps            = %d   (cap 2^33 = %d ; equal? %s)"%(J,2**33,J==2**33))
print("  dps/dpexp        = %.4f  (closed-form check, want ~1.0)"%(DP/DX))
print("  candidates       = %d"%C)
conf = 1-math.exp(-J/(2*2**30))
if C==0 and J==2**33 and abs(DP/DX-1)<0.02:
    print("  VERDICT: k0 NOT in [0, 2^60), at confidence %.1f%% (exponential model, mean 2*sqrt(L))"%(100*conf))
else:
    print("  VERDICT: NOT A CLEAN MISS -- do not cite. conf would have been %.1f%%"%(100*conf))
PY
  else
      echo "STATUS: no DONE line -- the run did not complete.  NOT a result."
      tail -1 head60b.err 2>/dev/null
  fi
} > R60_VERDICT.txt 2>&1
