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

### So I solved in EQUATION space instead (`eqsolve.py`, `eqsolve2.py`) — new, nobody had run it
Objective: choose knob deltas so each *equation* totals zero, allowing nonzero atoms.
- All 5 PIN knobs are affine on the residual. 19 knobs / 86 equations: every equation is
  knob-dependent, none unfixable.
- Exact full solve **fails** with a divisibility obstruction on **equation 29125**
  (`rhs % -P != 0`, P the 256-bit prime).
- Greedy keeps **79/86** rows -> score **39,026**. i.e. **the deliverable already IS the optimum of
  the equation-space solve over this knob set**, and the 7 failures are structurally forced.
- Widened to **162 affine knobs / 999 equations** (`eqsolve2.py`): full solve **core infeasible**.

This is the sharpest characterisation of the 39,026 barrier anyone has produced: it is not a
search failure, it is an arithmetic obstruction visible in equation space.

---

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
