# RESUME_M — agent M. Round 2: fix the representation, then price root-firing configurations.

## 0. Headline
1. **THE GATE IS PASSED.** E's engine could not represent the 39,026 deliverable. The defect is
   **structural, not numeric**, I found the exact cause, and the fix costs 5 variables.
   `engine2.py` now reproduces the deliverable **exactly** (zero vars differ), scores
   **39,026/39,033** with the same 7 failing equations, and keeps all **8 nonzero atoms** intact.
2. **I RETRACT part of my own §9 from last round.** I claimed E's enumeration never fired the
   root gate. That is **wrong** — E's reps for the 178-block were `47`/`112`, which are A-side,
   so E's `mask(1,·,·)` rows *are* root-firing. Details in §3. The coordinator recorded my claim;
   it needs correcting in FLEET.md.
3. **I did not beat 39,026.** Every search below baseline; all numbers reported.

---

## 1. STEP 1 — the representation defect, diagnosed and fixed

### What it actually was
Not the quadratic-branch ambiguity I expected. **All 23 differing vars are LINEAR, single-root.**
E's `_bootstrap` assigns each derived var `u` a *definer atom* `i`, and `forward` then **solves
that atom to zero**. So **any atom used as a definer is identically satisfied in every state E's
engine can produce.** Five of the deliverable's eight nonzero atoms are definers:

| atom | defines | at deliverable | at E's forward |
|---|---|---|---|
| 23616 | x_7068  | nonzero | **0** |
| 23617 | x_28730 | nonzero | **0** |
| 36659 | x_29854 | nonzero | **0** |
| 36663 | x_31864 | nonzero | **0** |
| 36664 | x_642   | nonzero | **0** |

Atom 36663 is literally the expression `x_31864` — E defines that variable as "= 0".
The other 18 differing vars are pure downstream contamination (all 18 have `delivIsRoot=True`:
their defining atom, evaluated at the deliverable, returns the deliverable's own value).

### The fix
`engine2.py`: demote those 5 atoms from definer role, promoting their 5 vars to free inputs.
`|FREE| 8365 -> 8370`, `|SEQ| 30383 -> 30378`. Nothing else changed.

### Gate result (`gate.py`)
```
PIN (vars promoted to free): [642, 7068, 28730, 29854, 31864]
GATE 1 exact vector reproduction : vars differing from deliverable = 0
GATE 2 score                     : satisfied 39026/39033 (failing 7)
         failing [12231,12270,12350,14584,18673,22044,29125]  == checker.py's list
GATE 3 the 8 nonzero atoms       : {23616,23617,36659,36660,36661,36662,36663,36664}  MATCH
GATE PASSED
```

### engine2 independently validated against `checker.py` on NON-deliverable points
Not just the deliverable — otherwise the fix could be overfitted:
| point | engine2 predicts | `checker.py` says |
|---|---|---|
| deliverable + leaf 4287 ON | 39,000 | **39,000** |
| deliverable + leaf 17378 ON | 38,961 | **38,961** |

### Bonus capability the fix unlocks
The 5 promoted vars are **the knobs that drive 7 of the 8 residual atoms affinely**. In E's
engine they were derived and therefore unavailable to any solver. They are now free knobs.

---

## 2. STEP 2 — pricing root-firing configurations

Deliverable ON-leaves: **24601 (A-side, block 178)** and **2081 (B-side, block 21)** — one per
root slot, so the root gate fires. cfg0's `{1530,1603}` are both B-side.

### The deliverable is a strict LOCAL MAXIMUM in leaf space (`enum2.py`)
| move from the deliverable | best score |
|---|---|
| turn OFF 2081 | 38,872 |
| turn OFF 24601 | 38,909 |
| turn OFF both | 38,776 |
| turn ON any 1 of the other 254 leaves | **39,000** (leaf 4287); histogram 38,961–39,000 |
| turn ON any 2 (top-40 pairs, 780 configs) | **38,975** |
| simultaneous repair (`simsolve`) at the base | **39,008** |

So both adding and removing live leaves cost. **Nothing reached 39,026.**

### The key methodological finding: atom-count is the WRONG objective
`simsolve` drives every bad atom to zero and *loses* (39,026 -> 39,008). The deliverable keeps
8 atoms nonzero **because they cancel inside the equations**: 8 nonzero atoms produce only
**7 failing equations**. Minimising atoms is not minimising failures.

### So I solved in EQUATION space instead — new, nobody had run it
Objective: choose knob deltas so each *equation* totals zero, allowing nonzero atoms.

### !! WITHDRAWN: there is NO divisibility obstruction on equation 29125 !!
I reported one at check-in 8. **It was wrong, and I withdraw it.** The message
`rhs % -P != 0` came from **one elimination ordering inside a badly overdetermined window**
(86 rows against 19 knobs), not from any invariant. Measured properly (`eq29125.py`, `eqsub.py`):

- Single-row solvability is exact and window-independent: `sum_f coef_f d_f = -s0` is solvable
  over Z **iff gcd(coef_f) | s0**. For equation 29125 the gcd is **1**. It divides everything.
  **Row 29125 is individually solvable**, and when actually solved it *is* zeroed.
- Same for all seven: gcd 1 for six of them, 40490 for eq 22044 — **all seven divide, all seven
  are individually solvable.**
- `eqsub.py` then removed the window question entirely: solve each subset of the 7 failures, then
  **apply it, re-propagate, and measure the true score**, so collateral damage is counted by
  measurement rather than assumed. Result: **all 127 non-empty subsets are solvable, 0 infeasible**,
  including all 7 at once, and in every case the solver really does zero its targets.

### What the barrier actually is: a minimum-cost residual, not an obstruction
Every repair fixes its target and breaks strictly more elsewhere:

| fix | failures 7 -> | score |
|---|---|---|
| eq 12350 | 10 | 39023 |
| eq 18673 | 10 | 39023 |
| eq 12270 | 11 | 39022 |
| eq 12231 | 18 | 39015 |
| eq 22044 | 28 | 39005 |
| eq 14584 | 34 | 38999 |
| eq 29125 | 34 | 38999 |
| **all 7** | **44** | **38989** |

Best over all 127 subsets: **39,023** — below the 39,026 baseline. So 39,026 is a **strict local
optimum in equation space**, and equation 29125 is not obstructed but merely one of the two
*most expensive* rows to repair.

### Why the residual is so cheap: it is confined to 7 equations
The 8 nonzero atoms are seen by exactly these 7 equations, nested, all with zero constant part:

    eq 12270 / 12350 / 14584 : all 8 atoms      eq 12231 : 6 of them
    eq 18673 : 36660,36661,36662,36663          eq 22044 : 23616,23617,36664
    eq 29125 : 23617 only

Equation 29125 sees a **single** nonzero atom, 23617 = `x_28730 - x_17499*x_9413`. The knob that
moves it directly is **x_28730 — one of the 5 vars my fix freed**. Zeroing it is precisely what
E's orientation did by construction, and it costs 27 extra failures.

### Answers to the questions asked
- **Q1 what divides what:** nothing obstructs — the gcd condition holds in all 7 cases.
- **Q2 which knobs move it:** eq 29125 has 12 affine knobs; the direct one is `x_28730`.
- **Q3 property of eq 29125 or of the window:** **of the window.** The row alone is solvable.
- **Q4 is "core infeasible" at 162 knobs instance or widening:** **the widening.** That window was
  **999 rows against 162 knobs** — 6:1 overdetermined, so generically infeasible regardless of the
  instance. Now proven rather than argued: every one of the 127 subsets is feasible.
- **Q5 knob set:** the 5 freed definer vars `[642, 7068, 28730, 29854, 31864]` **are** in it and
  **all 5 are affine**, verified in every run (`eq29125.py` prints the check per knob set). The
  widest set tried was every free var in the cone of every atom of the target equations.

---

## 2b. x_7068 identified — the deliverable corrupts exactly FOUR handles
P reports the four corrupted handles `x642, x28730, x29854, x31864` are exactly the four of 3,707
for which `P` does not divide the value. Those are four of my five freed vars. **`x_7068` is
collateral, not a fifth corruption:**

| | the four handles | x_7068 |
|---|---|---|
| definer form | product (`h = a*b`), or bare (`x_31864`) | **linear combination** `x_2099 + 7376877*x_642` |
| magnitude | 723, 89, 724, 724 digits | **90** digits (x_2099 is 89) |

Atom 23616 is 730 digits and the `-7376877*x_642` term alone is 730 digits;
`atom23616 + 7376877*x_642 = x_7068 - x_2099` is only 89. And for atom 23616 to be *satisfied*,
`x_7068` would have to be **730** digits — it is 90. So **the deliverable left `x_7068` at its
natural value and let the atom carry `x_642`'s corruption.** `P` does not divide `x_7068` either,
but that is automatic for a linear combination containing a corrupted term; `x_7068` is not
product-defined, so it would not sit in a handle population built from product definers. P finding
exactly four failures is consistent with `x_7068` not being in that population at all.
**Answer: a definer that is not a handle.** (Verified in my frame; I did not read P's handle list.)

## 2c. Alternative placements priced — none beat 39,026
`engine3.py` generalises the engine to an arbitrary demotion set. Demoting an atom and seeding its
variable with its current value is **bit-identical**, so demotion is score-neutral and purely adds
freedom (validated per candidate). Candidates: the **10** atoms that are currently zero, sit in a
failing equation, and are definers. Demoting one frees a variable that can *cancel* the existing
nonzero contribution inside that shared equation rather than zeroing atoms.

All 11 placements (baseline + 10 single demotions), equation targets of size 1 and 2, solved,
applied, re-propagated, scored: **every one returns 39,026; none above.**

**Structural constraint found:** equations **12270** and **18673** have **zero demotable zero
atoms**, so 2 of the 7 failures cannot be touched by this move at all.

**SCOPE — state it plainly.** `eqsub` prices *repairs of the current placement*. `place.py` prices
*alternative placements in a local neighbourhood of it* (one extra demoted atom drawn from the
current failing equations). **Neither prices a genuinely different placement** — corrupting a
*different set of handles*. That space is not reachable by adding atoms to the current failing
equations, and the handle population is P's object, not enumerable from the residual side.
**That is where anything above 39,026 would have to live, and it remains untested.**

## 2d. The candidate-agnostic pricer — built, calibrated, and used (round 5)
`price.py` + `engine3.py`. Input: **handle variables only**; the collateral demotion is derived.
From `[642, 28730, 29854, 31864]` the closure returns exactly engine2's PIN and DEMOTE sets, and
`price_given` returns **39,026, the same 7 failures, 0 vars differing**.

**My first tuner failed calibration (39,008 on the deliverable's own site) and I did not report its
twelve candidate numbers as results.** Diagnosis refuted my own hypothesis: the deliverable fixes
18 baseline failures and breaks **zero**, so there is no trade; and the affine model predicts
exactly at a delta of 10^728 (12/12). The model was sound, my **solver** was wrong — a ~40-knob
set let the sparse solver pick degenerate solutions. Restricting knobs to the freed handles fixes
it: `TUNER CALIBRATION 39008 -> 39026, PASSED`.

### The structural result (the useful half)
| | rows_target | tuned |
|---|---|---|
| deliverable's site | **25 of 25** | **39,026** |
| all 11 other L sites | **0** | 39,008 |

> **A site can help only if its corrupted atoms appear in the equations that fail at the
> uncorrupted baseline.** Those 25 equations are: 2554, 5324, 6816, 8124, 8680, 9041, 9123, 9421,
> 11226, 12231, 12270, 12350, 14584, 15558, 18673, 21000, 22044, 22534, 22997, 28929, 29125,
> 29330, 32026, 35512, 38051.

L's incidence measure is not this quantity — its top 12 are all 0-incident. Anything with
`rows_target = 0` can be discarded without pricing. (Caveat: the 25 are relative to this baseline.)

**Throughput:** `price_given` 0.53 s (~6,700/hour); tuned 1 s (0-incidence) to 4 s (fully
incident) → **~900–2,700/hour single-core**, 4 cores. List size is not the constraint.

## 2e. The filter is baseline-independent, and the incident set is 0.28% (round 6)

**Caveat resolved.** The 25 equations were computed in E's defective orientation. Recomputed
independently against the deliverable's own baseline — start from the deliverable's vector in the
corrected engine and un-corrupt it in place (each freed var back to what its own definer atom
prescribes, iterated to a fixpoint since x_7068's definer references x_642), then re-propagate:

| | score | failures | bad atoms |
|---|---|---|---|
| A: E's orientation, full forward | 39,008 | 25 | 5 |
| B: deliverable un-corrupted in place | 39,008 | 25 | 5 |

**IDENTICAL — `|A ∩ B| = 25`, no differences either way.** The filter is a property of the
instance, not of the orientation. **L can use it.** (`baseline_sets.json`)

*Caveat that remains real:* both baselines share the deliverable's FREE INPUTS. The 25 hold at
this free-input configuration, verified across two orientations; a materially different
configuration is untested.

**How thin is the incident set?** From the equation side only (no site enumeration):
moving handle `u` changes every atom mentioning `u`, so `u` is incident iff `occ[u]` meets the
78 atoms appearing in the 25 equations.

> **32 incident handles out of 11,307 (0.28%)** — and **all four** of the deliverable's are in it.
> `[642, 1627, 1844, 1956, 2218, 2892, 4863, 6480, 7062, 7945, 9629, 10861, 11425, 15422, 16495,
> 21279, 21718, 23538, 23642, 23754, 23822, 26732, 28098, 28730, 29305, 29854, 30175, 31465,
> 31864, 33001, 35619, 37413]`

Cross-check, both directions: L's 11 sites have **0 of 44** handles in the pool (matching their
measured `rows_target = 0`); the deliverable's four are **4 of 4** (matching `rows_target = 25`).

**Consequence:** the space collapses from C(11307,4) = 6.81e14 to **C(32,4) = 35,960**, a
**1.9e10** reduction — **~20 core-hours, ~5 h on 4 cores: exhaustively priceable.**
Incidence does *not* make the deliverable unique (it is one of 35,960); and the criterion is
necessary, not sufficient, so 32 is an upper bound.

## 3. E's monotonicity: ARTIFACT — but not for the reason I gave. I correct myself first.

### 3a. My §9 mechanism was wrong
I wrote *"the root gate never fires at cfg0 or at any configuration E enumerated from it."*
**False.** E's representatives for the 178-channel are `47` and `112`, both A-side; with
`1530`/`1603` (B-side) already ON, every `mask(1,·,·)` row puts leaves on **both** root slots.
E enumerated **12 root-firing configurations** (38,872–38,987).

### 3b. But the conclusion holds, with the right mechanism
Every root-firing configuration E priced **also carried cfg0's two B-side leaves**, i.e. it had
**≥3 live leaves**. **E never priced a clean 2-leaf root-firing configuration** (one A-side + one
B-side and nothing else) — which is exactly the deliverable's pattern. Its root-firing rows lost
because they paid the 3-live-leaf penalty, not because root-firing is bad.

### 3c. Measured from a NEUTRAL base — the curve is UNIMODAL, not monotone (`rfenum2.py`)
All 2^8 block-subsets, one representative per block, exact re-propagation, 256 configs in 13 s:

| n live blocks | 0 | 1 | **2** | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| best raw | 38776 | 38791 | **38804** | 38787 | 38773 | 38761 | 38746 | 38718 | 38655 |

It **rises** from 0 to 2 and falls after. "Monotone decreasing in live channels" is an artifact of
measuring only *additions* from a base that already sat at the 2-leaf peak.

**And the peak is root-firing:**
| | best raw |
|---|---|
| root-firing (A-side + B-side) | **38,804** |
| 78-side-only | 38,791 |
| **at exactly 2 live blocks:** root-firing | **38,804** |
| at exactly 2 live blocks: 78-side-only | 38,774 |

Best configuration in the whole lattice: `(47, 490)` — **47 is A-side (block 178), 490 is in
block 21** — the *same slot pattern as the deliverable* (`24601` A-side, `2081` in block 21),
found independently from a neutral base. Root-firing wins by **30 points** at equal channel count.

### 3d. The other, larger limitation: the instrument
Separately, `simsolve` overwrites the tuned handle values, so it cannot express a tuned-handle
optimum and tops out near 39,005 (cfg0) / 39,008 (deliverable base) from *any* base. The raw
numbers above are ~38,800 for the same reason — handles are untuned there. So 39,026 is out of
reach of every simsolve-based enumeration regardless of base.

**Answer to the coordinator's question: ARTIFACT.** Root-firing configurations are strictly better
at equal live-channel count, the score curve is unimodal with its peak at the deliverable's slot
pattern, and E's monotonicity came from enumerating only additions to an already-peaked base.

---

## 4. SECONDARY — slot-structure oracle for F: NEGATIVE, do not use the finer version
`orefine.py`: 136 bases (one saturating each tree stage) on the corrected engine, signature =
delta of each base's *own* bad-atom set. Maximal common refinement = **32 blocks**
(sizes 142,36,16,9,8,8,4,3,2×6,1×18) vs the old **8**. Four times finer — **but it is not valid.**

Validity, checked two independent ways (`ofinal.py`, `xslot.py`):
| partition | crossings vs tree96 (laminar test) | vs F's decoded slot pairs |
|---|---|---|
| old 8-block (E engine, ROWS sig) | **0** | 3 consistent, **0 cut across** |
| new 32-block (corrected engine) | **10** | 3 consistent, **2 CUT ACROSS** |

The two confirmed disagreements are **stage 27596** (F decodes 10\|12, oracle says 20\|2) and
**stage 27257** (F decodes 8\|3, oracle says 9\|2). The extra resolution is partly spurious.
The only stages the oracle splits that are absent from `mux_wiring.json` are `10136`, `15431`,
`21279` — each a trivial 1\|1 split of a 2-leaf stage, i.e. **no new information**.

**Recommendation to F, weakened by agent P's audit (see §4a).** Use the 8-block partition
(`blocks8.json`) as the trustworthy artifact. The 32-block partition is **unvalidated, not
refuted** — both yardsticks I judged it against (`tree96.json` and `mux_wiring.json`) are, per P,
coarsenings of the real tree, so its 10 "crossings" and 2 "cut across" verdicts may be the
yardstick's error rather than the oracle's. Files: `oracle_for_F.json` (both partitions +
validity flags + caveats), `xslot_report.json`. I did **not** modify F's directory.

### 4a. Scope: which of my claims are stated against `tree96.json`'s numbering
Agent P reports the real tree is **383 law-blocks at depth 9**, not 96 stages at depth 6, i.e.
`tree96` is a **coarsening**. Reading my results against that:

**Unaffected** (they never used tree96's node set):
- The 8-block refinement `178 | 41 | 21 | 6 | 3 | 3 | 3 | 1` — a measurement over the 256 leaves.
- The **root split 178 | 78**, which P independently confirms.
- Every score, the gate, `engine2`, the equation-space obstruction, and all of §2/§3.
- The block-for-block agreement with `mux_wiring.json` at the root, 19538, 10649, 21408 — those
  compare leaf **sets**, which stay valid however the stages are numbered or subdivided.

**Stated against tree96's numbering, therefore scope-limited:**
- "Of tree96's 88 non-root stages, exactly 2 are ever split; 86 are never cut." This is a
  statement about the **coarse** object. The fine tree has stages my 339 configurations never
  separated, so it is a limit on scope, not a contradiction — and it sits alongside the 8-block
  resolution limit I already reported.
- The stage IDs throughout (19538, 10649, 21408, 27596, 27257, 30973, 24533, 10136, 15431, 21279)
  are tree96 labels; a later reader must not assume they are fine-tree nodes.
- Consequently §4's negative verdict on the 32-block partition is **downgraded to unvalidated**.

Note my first crossing test was itself buggy — it flagged a block *containing* a whole stage as a
violation, which is legitimate for a laminar family. Corrected criterion: a crossing is an overlap
that is neither containment nor reverse containment. The numbers above use the corrected test.

---

## 5. Files added this round (all in `agentM_work/`)
`engine2.py` (corrected engine) · `fast2.py` (incremental, validated == full forward) ·
`fscore.py` (exact equation scorer, **2331x** faster than `E2.eqfails`, verified identical) ·
`mcore2.py` · `chan2.py` (corrected channels + simsolve) · `gate.py` (**the gate**) ·
`diag1.py`/`diag2.py` (the diagnosis) · `enum2.py` · `rfenum2.py` · `eqsolve.py`/`eqsolve2.py` ·
`orefine.py` · `ofinal.py` · `xslot.py`
Data: `blocks_corrected.json`, `oracle_for_F.json`, `xslot_report.json`, `enum2_raw.pkl`,
`enum2_pairs.pkl`, `orefine_sigs.pkl`, `diag1.pkl`, `val_4287.json`, `val_17378.json`.
Logs: `orefine.log`, `rfenum2.log`, `eqsolve2.log`, `pairs2.log`.

Round-1 findings (channel partition == tree, set-for-set) stand unchanged except the §9
retraction above. No other agent's directory was touched; no git commands were run.

## 6. Highest-value next step
Do **not** spend more cores on simsolve-based enumeration — it provably cannot express a
tuned-handle optimum. The live question is equation 29125's divisibility obstruction: find a knob
set in which that row becomes solvable, or prove the obstruction is invariant. `eqsolve2.py` is
the harness for it; it currently reports "core infeasible" at 162 knobs / 999 equations.

---

# ROUND 14 (post-restart #2) — LIVE STATE. Read this first.

## R14.0 The environment resets `*.pkl` globally. Rebuild before measuring ANYTHING.
`solve_lab/.gitignore` (and the repo root's) carry `*.pkl`, so every restart wipes the whole
parse/orientation cache chain, mine included. The rebuild is three commands and ~90 s:

    cd solve_lab/agentM_work
    python3 -u _rb_parse3.py   # copy of agentE_work/parse3.py -> model3.pkl  (atoms 40,727)
    python3 -u _rb_dag.py      # copy of agentE_work/dag.py    -> dag.pkl     (free 8,365, seq 30,383)
    python3 -u calib_r14.py    # imports shim -> harness_m.py; builds orient.pkl; runs the gates

`shim.py` must be imported BEFORE `ieng`/`price`, so `harness` resolves to `harness_m.py`
(E's harness with the three pkl paths repointed into agentM_work). Everything I run does this.
Do NOT write into agentE_work.

**Careful with backgrounding**: `cd X && cmd & echo $! > p.pid` backgrounds the whole `cd X && cmd`
chain, so the `.pid` lands in the ORIGINAL cwd. I put four pid files in the repo root that way and
moved them back. Use `nohup bash -c 'cd X && cmd' & echo $! > /abs/path/p.pid`.

## R14.1 Gate G1 re-passed after the rebuild (`calib_r14.log`)
    model: atoms 40727  eqs 39033  free 8365  seq 30383      (all pre-restart values)
    baseline: 39008, 25 failing, 5 bad atoms
    G1  witness {642,28730,29854,31864} -> 39026, fails exactly the 7,
        8 bad atoms, vars differing 0 of 38,748                 PASSED
    G2  39026 / 39000 / 38961 on the three CLI-agreeing points  PASSED
    G3  T's 12 cofactors zeroed -> 39021, 12 failing, list ==   PASSED
    G4  incremental == full engine3 (0 vars differing)          PASSED
    G5  tune() 39008 -> 39026 at nprobe=10                      PASSED
    G5b tune() 39008 -> 39026 at nprobe=80                      PASSED
    G6  0.006 s/site on general 4-subsets (box is uncontended now)
`checker.py solve_lab/best/new_instance_partial_39026.json` -> 39026/39033, failing
[12231,12270,12350,14584,18673,22044,29125]. Deliverable verified independently.

## R14.2 What was already done between LOG_M §81 and restart #2 (logs only, never logged)
    enumsub2.py      resumable rewrite: real resume, per-size distributions, errors COUNTED
    2^12             COMPLETE 4,096/4,096, above 39,026 = 0
    2^16 @ p10/30    COMPLETE 65,536/65,536, above 39,026 = 0, best 39,026 at the witness
    verifysup16      all 114 subsets at 39,026 over 2^16 contain the witness; 0 do not
    granul           p10/30 -> p80/180 moved 773 of 1,784 sampled subsets UP (max +10), 0 above
    2^16 @ p80/180   COMPLETE through |W| = 8, killed at |W| = 9
    2^18 @ p10/30    COMPLETE through |W| = 11, killed inside |W| = 12 (234k/262k)

## R14.3 Granularity: axis 1 is SATURATED, axis 2 is live
`tune()` probes `sols[j]`; `len(sols) <= |FAILS_UNC| = 25`, so at `nprobe = 80` the index set
already covers every solution. **`nprobe > 80` is a no-op — do not spend cores on p400.**
The live axis is the greedy ROW ORDER. `gran2.py`, 1,193 subsets x 9 orders: **1,001 moved up,
max +12, 0 above 39,026, witness optimum attained at the identity order.** Every per-subset score
in every distribution I have published is therefore a LOWER BOUND; the maximum is not.

## R14.4 The axis nothing has explored: the cofactors are not knobs
All 12 cofactors are FREE INPUTS, hence never in a closure, hence never in `tune`'s knob set.
The enumeration varies which handle relations break while holding the cofactors at the
deliverable's values. `enumcof.py` widens the knob set to closure(W) u cofactors (4 or 12).

## R14.5 Scripts added this round (all in agentM_work)
    _rb_parse3.py _rb_dag.py   rebuild the pkl chain after a restart
    calib_r14.py               the gates, shim-first          -> calib_r14.log
    gran2.py                   row-order granularity          -> r14_gran2_16.log, gran2_16.json
    enumsub3.py                enumeration, max over row orders (resumable)
    enumcof.py                 enumeration with cofactor knobs (resumable)
    xcheck14.py + r14_runchecker.sh   materialise 12 spread subsets and run checker.py on each

## R14.6 Jobs in flight (resume points if a third restart lands)
Every enumeration checkpoints to a `.pkl` every 2,000 subsets and replays from the stored index,
so a restart costs only the rebuild — but the `.pkl` itself is wiped by the restart, so in
practice a restart means re-running from 0. Budget accordingly: 2^16 is ~20-40 min under fleet
load, 2^18 is ~4x that.

    r14_enum16_d.log     2^16 @ p10/30  -- I KILLED it at index 40,000 with |W|=0..8 COMPLETE
                         (its full result is on record from the pre-restart run, and this run
                         reproduced all nine distributions entry-for-entry).  ckpt enumsub16.pkl
    r14_enum16_p80.log   2^16 @ p80/180 -- the headline run.  ckpt enumsub16_p80.pkl
    r14_enumcof16.log    2^16 @ p80/180 with the 4 cofactor knobs.  ckpt enumcof16_c4_p80.pkl
    r14_enum18_p80.log   2^18 @ p80/180 -- chained by r14_chain18.sh to start only once the
                         2^16 p80 run prints its final BEST.  ckpt enumsub18_p80.pkl
    r14_checker.log      12 materialised subsets, each scored by checker.py from outside my parse

If anything ever exceeds 39,026 the enumerators write the full assignment to
`M_sub*/M_cof*/M_ord*_<score>_<handles>.json` themselves and print `*** ABOVE 39026 ***`.
Grep the logs for `ABOVE` before trusting a summary.

## R14.7 RULE learned this round — the pricer is not monotone in the knob set
Adding knobs can LOWER `ieng.tune`'s score. Measured: at |W| = 5 over 2^16, 12 subsets reach
39,026 with handle knobs only but only 1 does with the 4 cofactor knobs added. Mechanism: more
columns make more rows individually solvable, the greedy therefore KEEPS a larger row set, and
the larger system's solution scores worse. `tune()` maximises over probes inside one greedy
chain, never over chains.

> **Take the maximum over knob sets and over granularities, never within one instrument.**
> A widened knob set is a different instrument, not a refinement.

Corollary: the cofactor run does NOT supersede the handle-only run, and neither supersedes the
row-order-varied run. All three are needed and the reported best is the max of the three.
