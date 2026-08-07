#!/bin/sh
cd "$(dirname "$0")"
while pgrep -x xmitm > /dev/null; do sleep 30; done
echo "ALL SIZE-5 SCAN PROCESSES EXITED at $(date -u +%H:%M:%S)"
grep "size=5" rep_real.txt
