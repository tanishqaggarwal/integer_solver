# CUSTODY.md — agent AI (custodian)

Custody opened **2026-08-07 21:25:43 UTC**. All times UTC.
Scope: keep the two orphaned jobs alive and honest, track disk, verify the restored
shared tables, audit the fleet for the two known false-record patterns.
**Read-only outside `solve_lab/agentAI_work/`.** No git commands run. Nothing deleted.

Artefacts in this directory:

| file | purpose |
|---|---|
| `verify_rot.py` | independent check of the rotational-sweep invariant (exit 3 on violation) |
| `verify_tbl.py` | independent check of the restored `tbl4s.bin` / `bm4.bin` |
| `audit_counts.py` | closed-form check of recorded scan counts; AA marker/engine cross-check |
| `sampler.sh` / `sampler.log` | 120 s resource + progress sampler (PID in `sampler.pid`) |
| `ALARMS.log` | written only when a threshold trips; absent = nothing tripped |

---

## 1. Job A — agent X's rotational sweep

**Status at 21:34: RUNNING. Not restarted — it never died.**

Identity established by PID, not by command-line matching:

```
/proc/30892/cmdline = /bin/sh ./rotall.sh
/proc/30892/cwd     -> /home/user/integer_solver/solve_lab/agentX_work
start (from /proc/30892/stat field 22) = 21:15:59
kill -0 30892 -> alive at 21:25:43, 21:27, 21:31, 21:32, 21:34
```

Measured rate: **8 rotations in 1086 s = 135 s/rotation** (X estimated ~130 s).
Remaining 120 rotations → **ETA ≈ 4.5 h from 21:34, i.e. ≈ 02:05 UTC.**

### 1.1 The correctness invariant — checked, holds

`verify_rot.py` recomputes, for every rotation with data, the sum of the six shard
candidate counts and compares with `C(128,5) = 264,566,400`. It also checks each shard
against **its own** closed form — the number of 5-subsets of `[0,128)` whose minimum
lies in `[lo,hi)` is `sum_{m=lo}^{hi-1} C(127-m,4)` — which is a stronger test than the
total alone, because two compensating shard errors would pass the sum but fail this.

Result at 21:32, rotations 0–7 (every rotation completed so far):

```
rot 0..7: shards=6  sum=264566400  delta=+0  per-shard-closed-form=all match  OK
rotations marked DONE AND invariant-verified : 8  [0,1,2,3,4,5,6,7]
rotations marked DONE but NOT verified       : 0
rotations with data but invariant VIOLATED   : 0
HIT lines: 0        total of all zero= fields: 0
```

No short sum, so no rotation needs re-running. `HIT` count is 0 — **no solution found
so far.** Re-run `python3 verify_rot.py` after any further rotations; exit 3 = violation.

### 1.2 The partial result, stated as X stated it

What `r` of 128 completed rotations excludes is

> `{ S : |S| = 10, ∃ j completed with |S ∩ A_j| = 5 }`

and **nothing more**. It is **not** "`r/128` of `|S| = 10` exhausted".

**`|S| ≤ 10` is claimable only when all 128 rotations finish.** A set such as `{0,…,9}`
is covered by exactly one rotation (`j = 5`), so no proper subset of the rotations
exhausts the family. If the sweep dies partway, the honest statement is the conditional
one above. As of 21:34 the correct statement is:

> every `|S| = 10` ON-set having a balanced 5|5 split at one of rotations 0–7 is
> excluded; `|S| ≤ 10` is **not** yet exhausted.

`agentX_work/rotprog.py` already prints exactly this and refuses the stronger claim
while `len(done) < 128`. It is the right tool; keep using it.

### 1.3 Robustness defects in the sweep (reported, not touched)

1. **`agentX_work/rotall.pid` holds the wrong PID.** It contains `30889`, which is the
   *launching bash* (`/bin/bash -c … setsid nohup ./rotall.sh … & echo $! > rotall.pid`);
   the real `rotall.sh` is `30892`. `$!` captured the `setsid` wrapper. Consequence:
   agent X left a watcher (PID 9890) looping on
   `until … || ! kill -0 $(cat rotall.pid)`, so it is watching a process that is **not**
   the sweep. It will not notice the sweep dying. (I hit the identical `$!`/`setsid`
   bug in my own sampler and fixed it by having the script write its own `$$`.)
2. **`rotone.sh:19,44` — `while pgrep -x xrot …`.** Name-based process identification.
   Any stray `xrot` anywhere on the box stalls the sweep silently and for ever.
3. **`rotall.sh:4`** only skips rotations marked `DONE`; a rotation marked
   `INCOMPLETE (n/6)` is *not* retried within the same pass. It is retried only if
   `rotall.sh` is restarted. That is the intended restartability, but it means an
   incomplete rotation sits unfixed until someone restarts.
4. **`rotall.sh:6-7`** decides `DONE` from `grep -c "rot=$J "` — a *line count*, with no
   exit-code test and no check of the candidate sum. It is materially better than a bare
   marker (it counts engine-written lines, and a dying shard writes no `DONE` line at
   all), but it would accept six lines with short counts. That gap is what
   `verify_rot.py` closes; the counts do in fact verify.

---

## 2. Job B — agent AB's `d_reg(4)` rank computation

**Status at 21:34: RUNNING, no read-off yet. Not restarted (and must not be).**

```
/proc/6881/cmdline = python3 ab_dreg3.py 2 4
/proc/6881/cwd     -> /home/user/integer_solver/solve_lab/agentAB_work
elapsed 3067 s (51 min) at 21:34, state R, ~51 % CPU
RSS 3046 → 3382 → 2972 MB across 21:25–21:34 (oscillating, not a monotonic leak)
```

`dreg3.log` last grew at **20:44:40** and ends mid-step:

```
n=4 : 18 vars (4 boolean), 23 generators
   d=2 :    30 rows x    39 support-cols   rank(M)=29   rank(M+targets)=32   not yet
   d=3 :   398 rows x   456 support-cols   rank(M)=364  rank(M+targets)=367  not yet
   d=4 :  3365 rows x  3247 support-cols   rank(M)=2838 rank(M+targets)=2841 not yet
   d=5 : 21057 rows x 17091 support-cols = 3.60e+08 cells (build 0s)
                                     <-- no result line yet
```

The quiet log is expected, not a hang: the line for a `d` is written only when its rank
finishes. The comparable `n=3, d=5` step (`6494 × 4679`, 3.04e7 cells) took 81 s; this
one is ~12× the cells with ~3.6× the rank, so a much longer elimination, further slowed
by the box being at load ~17 on 4 cores. The process is in state `R` and burning CPU.

**Already established by this run** (visible in the log, no interpretation added):
`n=2` solving degree **4**; `n=3` solving degree **5**. So the sequence is `4, 5, ?`.

### 2.1 AB's read-off, to be applied exactly as written

AB wrote this in advance so it could not be rationalised afterwards:

- **`ALL SELECTORS PINNED`** ⇒ `d_reg(4) = 5`, sequence **4, 5, 5** — growth
  **sublinear**, and **AB's §9.12 verdict re-opens.** AB identified this as the only way
  its own conclusion could be wrong. **If this is the outcome it is the most important
  thing on the fleet and must be reported to the coordinator immediately, ahead of
  finishing anything else.**
- **`not yet`** ⇒ `d_reg(4) ≥ 6`, sequence **4, 5, ≥6** — strictly increasing at every
  measured step.

`d = 6` is over this box's cap either way, so `≥ 6` is the strongest positive statement
obtainable here.

### 2.2 Standing instruction

**Do NOT restart this job if it dies.** It is expensive and the coordinator may not want
it re-run. If it dies, report to the coordinator with the tail of `dreg3.log`.

---

## 3. Independent verification of the restored shared tables

X rebuilt `tbl4s.bin` and `bm4.bin` after the deletion incident and reported them
bit-for-bit identical to the pre-deletion values. **Checked independently at 21:28–21:29
— every claim confirmed.**

```
tbl4s.bin  1,420,712,448 bytes   md5 3065a6f304bad45561d051f518b604a6   MATCH
bm4.bin      536,870,912 bytes   md5 f3e458ee2564f18eb20c25492390fa8b   MATCH
key count  1420712448 / 8 = 177,589,056                                 MATCH
first two  [208528404822, 231390034609]                                 MATCH
last two   [18446743699321287810, 18446743880247473500]                 MATCH
```

Two structural checks X did not claim, added as corroboration:

- file size is an exact multiple of 8 and `size/8` equals the claimed key count, so the
  count is not an independent assertion but forced by the file itself;
- keys are ascending on a stride sample (200,214 of 177,589,056, step 887) — a truncated
  or interleaved rebuild would show an order break;
- `bm4.bin` is exactly `1<<29` bytes, which is the size `xrot.c:143` hard-checks with
  `FATAL … exit(2)` before scanning.

**X's restoration is correct.** `verify_tbl.py` re-runs all of it; exit 3 on failure.

---

## 4. Disk and load

**Not an emergency. No deletion recommended at this time.**

| time | free MB | agentAA_work MB | agentX_work MB | load |
|---|---|---|---|---|
| 21:28:43 | 9,801 | 11,275 | 5,853 | 15.6 |
| 21:30:43 | 11,208 | 11,275 | 4,444 | 20.6 |
| 21:32:44 | 11,617 | 11,275 | 4,041 | 17.3 |
| 21:34:35 | 10,490 | 11,275 | 5,159 | 16.8 |

Free space **oscillates, it does not trend down.** The sawtooth is job A: each rotation
writes six `rt_*.bin` build chunks (~1.4 GB total), merges them into `xrot_tbl.bin`
(943 MB) plus `xrot_bm.bin` (512 MB), then deletes all of them (`rotone.sh:30,45`).
`agentX_work` therefore breathes between ~4.0 and ~5.9 GB with a period of one rotation.
**Trough observed so far: 9,801 MB free.** The sweep is in disk steady state, so 120 more
rotations do not imply 120 more GB.

The two standing consumers:

- `agentAA_work` **11,275 MB, flat** across the whole window — it is the largest single
  directory on the box but it is *not growing*. Mostly `tbl/v0..v7/t4` shards plus
  `tbl/bm4.bin`.
- `agentX_work` 4.0–5.9 GB, of which `tbl4s.bin` (1.4 GB) + `bm4.bin` (512 MB) +
  `sbm.bin` (512 MB) are persistent.

**Warning thresholds are armed, not tripped.** `sampler.sh` writes to `ALARMS.log` if
free space goes under 6 GB, if `MemAvailable` goes under 1 GB, if either job dies, or if
a `HIT` ever appears. `ALARMS.log` does not exist as of 21:34, which is the all-clear.

Memory is the tighter resource, not disk: 16,075 MB total, no swap. `MemAvailable` was
11,198 MB at 21:34 and job B alone holds ~3 GB. An OOM kill would take job B — the
expensive, non-restartable one — so `MemAvailable` is sampled every 120 s.

**Recommendation only, per rule 3 — I have deleted nothing and propose nothing be
deleted.** `tbl*.bin` and `bm*.bin` are fleet property; `agentAA_work/tbl/bm4.bin` and
`agentX_work/{tbl4s,bm4}.bin` are on the live data paths of AA, Y and the orbit scans.
If space ever does become tight, the only candidates I would even raise are X's
per-rotation scratch (already self-cleaning) and `DEAD_spart*.log`, and that is the
coordinator's decision, not mine.

---

## 5. Audit for the two known false-record patterns

Read-only sweep of every agent directory. **No fabricated record was found.** Every
recorded number I could test against a closed form is correct:

- job A rotations 0–7: six shards each, sum exactly `C(128,5)`, every shard matching its
  own closed form (§1.1);
- agent Y `rep_comp.txt`: `n = 32,640 / 2,763,520 / 174,792,640 / 8,809,549,056` for
  sizes 2/3/4/5, equal to `C(256,2..5)` exactly;
- agent AA: in all 64 `runs*/…txt` files the shell-written `SHARD` marker count never
  exceeds the engine-written `DONE` line count (24 = 24 on every finished tag, 20 = 20 on
  the in-flight one), and within each tag+size all eight table-shards report identical
  `n`.

What I did find is **unsound mechanism** — code that would record a failure as a success
if one occurred. Ranked by exposure. **I have not edited any of these.**

### Pattern (i) — marker written without testing an exit code

| # | file:line | detail |
|---|---|---|
| **1** | `agentAA_work/aa_shard.sh:11-12` | **Highest severity — this script is running now (PID 13873).** Engine output is fully masked (`./aa_signed scan … >/dev/null 2>&1`) and `echo "SHARD$s" >> "$out"` is unconditional. Worse, the resume guard at **line 10** (`grep -q "sz=$b .*SHARD$s\$" … && continue`) keys on *that same shell-written marker*, so a shard that segfaulted would be recorded as done **and permanently skipped on restart**. This is agent Y's failure mode reproduced exactly. No damage yet — markers and engine `DONE` lines currently agree everywhere. |
| 2 | `agentY_work/yrun.sh:12-17` | `./ymitm scan … ; echo "finished size $SZ …" >> yrun.status` — unconditional, for sizes 2,3,4 and 5, plus `ALLDONE`. The sibling `yorbit_run.sh:20-28` was fixed (`RC=$?`, `[ $RC -eq 0 ] &&`, plus a table-exists precheck) — `yrun.sh` is the unfixed twin. **The records happen to be true** (verified vs `C(256,b)` above); the mechanism is not. Note `yrun.status` lines 1–2 are both stamped `20:03:52`, the same-second signature of the Y incident — here it is **genuine**: the engine logs show size 2 at 0.0 s and size 3 at 0.3 s. |
| 3 | `agentU_work/u22_run.sh:9-10` and `u25_run.sh:6-7` | `wait` then unconditional `ALL_SHARDS_DONE` / `ALL_DONE`; `wait`'s status never inspected, so any of the four shards dying is invisible. |
| 4 | `s10/bn_key4_drive.sh:9-12` | `timeout 900 … & … wait; echo "round $start complete"` … `echo ALLDONE`. Neither `wait` nor `timeout` is checked — a `timeout` kill (exit 124) reads as a completed round. |
| 5 | `agentX_work/runAB.sh:3-8`, `runB.sh:7-10`, `scanphase.sh:9-17`, `runbsgs.sh:5-10` | Unconditional `=== … COMPLETE ===` / `DONE` echoes after unchecked engine invocations; `runbsgs.sh` additionally sends engine stderr to `/dev/null`. Historical (not running). `runAB.sh:6,8` is the sharpest case — it prints "`|S| <= 8 exhausted`" / "`|S| <= 9 exhausted`", i.e. a *mathematical claim*, with no test that the scan succeeded. |

### Pattern (ii) — output piped or masked, so failure reads as success

| # | file:line | detail |
|---|---|---|
| 6 | `agentAA_work/aa_run.sh:16-19`, `aa_run6.sh:16-19`, `aa_deep.sh:12-14` | `./aa_signed scan … 2>&1 \| tail -1` — the pipeline's exit status is `tail`'s, always 0 — followed by unconditional `OFFSET_DONE` / `DEEP_DONE`. Structurally identical to X's `gcc … \| head -3 && echo rebuilt`. |
| **7** | `agentT_work/t_rebuild.sh:5-16`, `t_rebuild2.sh:1-8`, `t_rebuild3.sh:4-17` | **Subtle and worth singling out: these look protected and are not.** All three open with `set -e`, but every build step is `python3 -u X.py \| tail -3`. A pipeline's status is the last command's, so `set -e` sees `tail`'s 0 and never fires. A crashing `buildall.py` / `calib2.py` still reaches `echo "=== REBUILD DONE"`. `set -o pipefail` is absent from all three. `t_rebuild3.sh:8-10` additionally uses `[ -f X.pkl ] \|\|` guards, so a stale `.pkl` from a failed earlier run is silently accepted as current. |
| 8 | `agentB_work/bloop.py:30` | `subprocess.run([… 'borient7.py' …], capture_output=True)` with the returncode discarded, then `pickle.load(open(W+'orient7.pkl','rb'))` on the very next line. If the child crashed, `orient7.pkl` is the **previous round's file** and the loop proceeds on stale data — X's stale-binary shape, in Python. |
| 9 | `agentN_work/poly68b.py:172-186` and `gb16.py:43-47` | Singular's returncode discarded. `TimeoutExpired` is handled explicitly (`'TIMEOUT'`), but a non-zero exit yields empty stdout, and `poly68b.py:186` writes that empty string into `runs/poly68b.json` as the `singular` field — a crash recorded as a result rather than as an error. |
| 10 | live command, agent AE (no file) | AE's build ran `gcc -O3 -march=native … -o aekang aekang.c -lpthread 2>&1 \| head -30 && echo BUILD_OK && ls -la aekang` — **verbatim the pattern that produced X's fabricated candidate count.** This instance happened to succeed: `aekang` is 21:25:47, newer than `aekang.c` at 21:25:36, and executing it prints its usage line. Reported because the pattern, not the outcome, is the hazard. |

### Not defects (checked and cleared)

`agentAE_work/ae_lib.py:82` records `rc=pr.returncode` and a stderr tail into its result
dict and prints them — correct. `agentX_work/xstest.py:40` passes `check=True`.
`agentY_work/ycheck.py:27` discards the returncode but validates the parsed output
against independently computed expectations, so a crash fails the check.
`agentR_work/bench.py:21` discards z3's returncode but records `'ERR'` on empty output —
degraded, not false.

---

## 6. Handover / what to do next

1. Re-run `python3 verify_rot.py` after each batch of rotations. Exit 3 = a rotation's
   shard sum is short; that rotation must be **re-run, not marked done**.
2. Watch `ALARMS.log`. Its absence is the all-clear.
3. On job A dying: it is restartable and skips rotations already marked `DONE`
   (`rotall.sh:4`). Restart it, and record here that it died, when, and why, with
   evidence. Note the wrong PID in `rotall.pid` (§1.3.1) — use `30892`, or re-derive.
4. On job B dying: **do not restart.** Report to the coordinator with the tail of
   `dreg3.log`.
5. On job B printing `ALL SELECTORS PINNED`: report to the coordinator **immediately**,
   ahead of everything else — `d_reg(4) = 5`, growth sublinear, AB's §9.12 re-opens.
6. Do not claim `|S| ≤ 10` until `rotdone.txt` shows all 128 rotations `DONE`
   **and** `verify_rot.py` exits 0 on all 128.

### Custody log

- **21:25:43** custody opened. Both orphans confirmed alive by `kill -0` + `/proc/<pid>/cmdline`.
- **21:27** job A invariant verified for rotations 0–5.
- **21:28–21:29** restored tables independently verified (md5, key count, boundary keys, order, bitmap size).
- **21:29–21:32** fleet-wide read-only audit; findings in §5. No false record found.
- **21:32** job A invariant verified for rotations 0–7 (all completed rotations).
- **21:34** sampler restarted with memory tracking and alarm thresholds (PID 2490).
- No process killed except my own sampler (12126, replaced by 2490). Nothing deleted. No git command run.
