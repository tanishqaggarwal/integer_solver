# Quarantine notice — `ALARMS.log.UNRELIABLE`

**Raised by:** agent AI (custodian), against its own monitor.
**Quarantined:** 2026-08-08 00:10 UTC.

## What is in the quarantined file

40 lines, all of the form `ALARM dreg(6881) DEAD`, stamped every 120 s from
**22:31:30** to **00:06:44**. **Every one of them is false.**

## What actually happened

Job B — agent AB's `d_reg(4)` rank computation, PID 6881,
`python3 ab_dreg3.py 2 4` — **completed normally at 22:29:50 UTC**, wrote its
read-off, and exited. The terminal lines of `solve_lab/agentAB_work/dreg3.log`:

```
     d=5 :  21057 rows x  17091 support-cols = 3.60e+08 cells (build 0s)
        rank(M)=15947  rank(M+targets)=15950   not yet   (6305s)
     d=6 : 105720 rows x  72858 support-cols = 7.70e+09 cells (build 6s)
        OVER CAP (4.0e+08) -- cannot decide this degree on this box
     -> n=4 solving degree >=6   (6330s)

MEASURED SOLVING DEGREES: {2: '=4', 3: '=5', 4: '>=6'}
```

A `6305s` timing and a rank pair only a real elimination can produce. Nothing was
lost. This is the **`not yet`** branch of the read-off AB wrote in advance:
`rank(M) = 15947` and `rank(M+targets) = 15950` differ, so the selectors are *not*
pinned at `d = 5`.

## The defect (mine)

`sampler.sh` named the **outcome** from `kill -0 6881` alone:

```sh
kill -0 6881 2>/dev/null && B=alive || B=DEAD
[ "$B" = DEAD ] && echo "$TS ALARM dreg(6881) DEAD" >> ALARMS.log
```

`kill -0` returning false establishes exactly one thing: the process is gone. It is
**equally consistent with success and with failure.** Turning it into `DEAD` is a
claim about the outcome, made without reading the evidence — and the evidence was
sitting in `dreg3.log` the whole time.

> **Absence of a process is not evidence of failure.**
> **A monitor must read the artefact before it names the outcome.**

This is the same rule that kept agent Y's lying status file out of the record,
applied to the thing that does the watching. It is the fleet's own failure mode
landing on the monitor: a status marker written without checking the thing it claims.

A second, compounding defect: the alarm was **not latched**, so one wrong inference
was restated 40 times and looked like 40 pieces of evidence. It was one.

## Why the lines were kept

They are a true record of a real defect and are worth keeping as such. Deleting them
would erase the evidence of the failure. They are quarantined, not removed — the
header of `ALARMS.log.UNRELIABLE` states on its first line that every line in the
file is false.

## What was fixed

1. `ALARMS.log` now contains only alarms that survive scrutiny. As of this writing
   that is exactly one line — `23:37:47 ALARM memAvailable=4MB < 1024` — which is a
   **true** reading of a real near-OOM moment during the `d=6` build, and is retained.
2. `sampler.sh` now distinguishes `COMPLETED` from `DIED` by **reading the job's log
   for its terminal line** before naming either, and quotes that line in the alarm.
   `DIED` is emitted only when the terminal line is *absent*.
3. Alarms are **latched** — a terminal condition is reported once, not every cycle.
4. **Job B is no longer watched at all.** It is finished and its result is recorded.

## Status of Job B — not lost, and recorded

`d_reg(4) ≥ 6`. Sequence **4, 5, ≥ 6** — strictly increasing at every measured step.
`d = 6` is over this box's cap (`7.70e+09` cells vs a `4.0e+08` cap), so `≥ 6` is the
strongest positive statement obtainable here, exactly as AB anticipated.
**AB's §9.12 verdict is NOT re-opened.**
