# LOG_M — agent M. "Is E's channel partition the same object as F's tree?"

Answer: **YES, exactly, set-for-set.** Confirmed from two sides that never touch each other's data.

---

## 1. Baseline, re-verified by me first
```
python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json
[checker] satisfied 39026/39033  (7 failing)  eval=0.3s
failing [12231, 12270, 12350, 14584, 18673, 22044, 29125]
```
Atom level (E's model, `engine.badatoms` on the full deliverable vector): exactly **8 nonzero
atoms** `{23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664}`.
Deliverable ON-leaf set: `{24601, 2081}`. **I did not beat 39,026.**

`orient.pkl` regenerated in my directory is **byte-identical** to E's, so every measurement below is
in exactly E's model, not a lookalike.

---

## 2. The bug in the premise (and why the prediction looked untestable)
E: "16 booleans moving nothing at all". Reproduced `channels(triple8_seed)` exactly — same
178/41/21 sets as `chan_cfg0.json`, `exact same partition: True`. But the 16 are

> **14 genuinely inert + 2 already ON (`x_1530`, `x_1603`).**

`channels.py` line `if v0[f]==1: continue` skips ON bits, so they fell out of every class and were
counted with the inert ones. And the two ON bits sit exactly where saturation predicts:

| ON bit | tree node it sits in | size | the other leaves of that node |
|---|---|---|---|
| `x_1603` | **19538** | 9 | all 8 measured inert |
| `x_1530` | **10649** | 4 | all 3 measured inert |

So "14 inert" was never a property of those bits — it was E's own saturation law
(one live leaf saturates its subtree) observed and not recognised.

---

## 3. The test E asked for. RESULT: they become live.
Knob set: the **256 boolean leaves** returned by `mcore.bools()` — the booleans of the cone of
ROWS `[7389, 10187, 20212, 20215, 28647]`. Each probed singly, 0 -> 1.

| selector configuration | what the 16 do |
|---|---|
| cfg0 (`triple8_seed`, 1530=1603=1) | 14 inert + 2 ON |
| cfg0 with **1530 = 1603 = 0** | **ch9 = node 19538 exactly; ch4 = node 10649 exactly; 3 more join a ch24** |
| deliverable cfg (`{24601, 2081}` ON) | **a single live channel of exactly 16** |

All 16 become live. **Prediction CONFIRMED.**

---

## 4. The partition IS the tree — measured over 339 configurations
`sweep1.py` (42 configs) and `refine2.py` (297 configs: all-off, each of the 256 leaves alone ON,
40 random pairs). For each configuration: measure the channel partition, then intersect every one
of `tree96.json`'s 88 non-root stages against it.

- At cfg0: **every one of the 96 stages lies entirely inside one class. Zero split.**
- Across all 339 configurations, the only stages ever split are **19538 (always 6|3)** and
  **10649 (always 3|1)** — and those are the stages' own two input slots, i.e. resolution
  *below* `tree96`'s granularity, not a contradiction.
- **Zero arbitrary crossings anywhere.** An arbitrary cut would move between configurations;
  these two never do.

(ON leaves cannot be flipped, so a class boundary at an ON leaf is a measurement artifact and is
excluded from the crossing count. Stated because including them inflates the count to 17/42.)

### The maximal common refinement: 8 blocks
Over all 297 configurations, two leaves stay together iff their signatures agree in every
configuration where *neither* is ON (`refine2.py`; the naive version is contaminated because each
leaf is uniquely ON in exactly one config, which splits all 256).

> **178 | 41 | 21 | 6 | 3 | 3 | 3 | 1**  (`blocks8.json`)

That is the finest partition E's residual cluster can *ever* resolve. 41+21+6+3+3+3+1 = **78**.

---

## 5. Independent cross-check from the circuit side — exact, set-for-set
Computed slot leaf supports from **F's `mux_wiring.json` only** (selector variables of each slot,
their defining atom's cone, intersected with the 256 leaves) — no use of my measurement:

| stage | circuit inA / inB leaf support | residual-side block |
|---|---|---|
| **15298** (root) | 178 / 78 | inA **== block 178, set-for-set**; inB == union of the other 7 blocks |
| **19538** | 6 / 3 | **== blocks 6 and 3, set-for-set** |
| **10649** | 3 / 1 | **== blocks 3 and 1, set-for-set** |
| **21408** | 10 / 11 | at the deliverable's config: live class of 11 **== inB exactly**; the ON leaf 2081 is in inA; the inert-inside-21408 **== inA \ {2081} exactly** |

Two machineries, no shared data, identical sets. **E's channel model and F's tree are one object.**

---

## 6. Mechanism, stated with knob set and configuration
Knob set: the 256 boolean leaves of `mcore.bools()`. Configuration: any (tested at cfg0,
all-off, 256 single-leaf, 40 random pairs, the deliverable's).

> **The channel partition at configuration C is the leaf-support partition induced by the tree's
> slot structure, cut at the deepest gate that C saturates.**
> Switching one leaf ON inside a subtree S saturates S and splits S's class into S's two slot
> supports; every other leaf of S then measures as inert.

Corollary (the residual form of "a gate fires only when both inputs are live): with only a-side
leaves ON, the whole 178-block is inert in all 5 rows; it becomes a live class the moment any
b-side leaf is ON. Measured at `alloff` (178 signature all-zero) vs `one_2081` (178 live).

The signature vectors at all-off nest strictly by support, one row added per level:
`41`: row 20212 only -> `24`: +20215,+28647 -> `9`: +10187 -> `4`: +7389.
Two blocks at the same level share the added row's value exactly (6 and 3 share the 10187 value;
3 and 1 share a different one). That is E's residual congruences resolved **per subtree** — the
payoff the task named.

---

## 7. Boundaries of the oracle (measured, so nobody over-claims it)
- **Resolution limit: 8 blocks.** The 5-row residual cluster cannot see inside the 178, inside the
  41, or inside 21408 beyond one slot split. It resolves the b-side spine and "a-side live / not".
- **Widening the signature does not help.** Using the full bad-atom delta *support* instead of the
  5-row projection splits all 256 leaves — each leaf has private pin atoms (`fullsig.py`).
- **The stage-wire oracle fails in E's engine.** Recording which of the 96 stages' wires change
  (`ancestors.py`): at all-off, 219 of 256 leaves change **no** stage wire; at cfg0, 229 change
  exactly one (the root). E's forward propagation does not realize the circuit's *intermediate*
  stage values — it is a propagation-with-defect, valid at the cluster residual only. So the
  residual side can never, by itself, decode a-side internal wiring. F's decode is not replaceable.

---

## 8. What I tried for score, and what it cost
- **E's `simsolve` from the deliverable's configuration: 39,008.** Its knob discovery is anchored
  on cfg0's cluster; at the deliverable the residual atoms under E's orientation are
  `{23618, 34120, 36660, 36661, 36662}` — a different cluster.
- **E's `forward` from the deliverable's free inputs: 39,008, not 39,026.** 23 derived variables
  differ; E's orientation zeroes the 8 atoms the deliverable deliberately leaves nonzero and pushes
  the defect further out. **E's engine cannot represent the deliverable's configuration.** This is
  the single most important negative result here.
- `scan_single.py` (each leaf alone ON, through `simsolve`): 38,842. Abandoned.
- `pairscan.py` (deliverable's 30 handle settings kept, leaf pair swapped over 14x78 pairs):
  best **38,944** at `(21266, 2081)`; the deliverable's own pair gives 39,008. The handles are
  tuned to the leaf pair, so swapping leaves without re-solving handles is a dead instrument.

---

## 9. The consequence nobody had noticed
E's whole channel enumeration (`chanenum.py`, "the empty set wins at 39,005", monotone in live
channels) is anchored at cfg0 — **whose ON-set `{1530, 1603}` lies entirely in the 78-side**
(19538 and 10649 are both under the root's inB slot). The root gate therefore **never fires** at
cfg0 or at any configuration E enumerated from it. The deliverable's ON-set has one leaf in each
root slot and scores 21 higher. **E's monotonicity result does not price the configurations that
matter**, and the "more live channels is always worse" conclusion is a statement about the
78-side-only regime, not about the instance.

---

## 10. Files
`mcore.py` (channel measurement at any seed, full 5-coord signature) ·
`xcompare.py` (class -> maximal tree nodes; crossing counter) ·
`sweep1.py` -> `sweep1.pkl` · `refine.py`/`refine2.py` -> `refine2.pkl`, `blocks2.pkl`,
`blocks8.json` · `fullsig.py` · `ancestors.py` -> `anc_alloff.pkl` · `mchan.py` (E's channels.py
with absolute paths, read-only copy) · `pairscan.py` -> `pairscan.pkl` · `scan_single.py` ·
`cfg0_M.json`, `deliv_seed.json`, `orient.pkl` (identical to E's).
E's and F's directories were not modified. No git commands were run.

---
---

# LOG_M ROUND 2 — fix the representation, then price root-firing configurations

## 11. The gate (step 1): E's engine could not represent the deliverable. Now it can.

### 11.1 Diagnosis — structural, not numeric
I expected a quadratic two-root branch ambiguity (`_solvevar` sets `v[u]=0` when
`len(rts)!=1`). **That hypothesis was wrong.** `diag1.py`: all 23 differing vars are
`kind=lin`, `nroots=1`. Zero are quad-defined.

The real mechanism (`diag2.py`). `harness._bootstrap` gives each derived var `u` a *definer
atom* `i`; `forward` then solves that atom **to zero** to get `v[u]`. Therefore **every definer
atom is identically zero in every state E's engine can reach.** Five of the deliverable's eight
nonzero atoms are definers:

    atom 23616 -> defines x_7068     atom 36663 -> defines x_31864
    atom 23617 -> defines x_28730    atom 36664 -> defines x_642
    atom 36659 -> defines x_29854

    value at the deliverable          value at E's forward
      23616  -155122693640347919...     0
      23617  1193251489541137844...     0
      36659  2793021083152746734...     0
      36663  -265159771897697970...     0
      36664  4033004069730628753...     0

`atoms[36663]` is the expression `x_31864` — a one-variable atom, so E's orientation literally
defines that variable as "= 0" while the deliverable sets it to -2651597718976...

The remaining 18 differing vars are downstream contamination: all 18 satisfy `delivIsRoot=True`
(their defining atom, evaluated at the deliverable's own vector, returns the deliverable's
value), so fixing the 5 roots fixes them by propagation. It did.

### 11.2 The fix and the gate
`engine2.py` demotes those 5 atoms from definer role -> their 5 vars become free inputs.
`|FREE| 8365->8370`, `|SEQ| 30383->30378`. Seed grows from 32 to 37 nonzero free inputs.

    PIN = [642, 7068, 28730, 29854, 31864]
    GATE 1  vars differing from deliverable : 0
    GATE 2  satisfied 39026/39033, failing [12231,12270,12350,14584,18673,22044,29125]
    GATE 3  nonzero atoms == {23616,23617,36659,36660,36661,36662,36663,36664}
    GATE PASSED

### 11.3 Guarding against an overfitted fix
A fix validated only on the deliverable proves nothing. Two independent non-deliverable points,
predicted by `engine2` then checked with the real `checker.py`:

    deliverable + leaf 4287 ON  : engine2 39000  ->  checker.py  satisfied 39000/39033
    deliverable + leaf 17378 ON : engine2 38961  ->  checker.py  satisfied 38961/39033

`fast2.resid_delta` also verified equal to a full `engine2.forward` + `badatoms` recomputation.

### 11.4 What the fix unlocks
The 5 promoted vars drive 7 of the 8 residual atoms affinely (`av[36663] = x_31864` exactly,
`av[36664] = x_642 - x_28599*x_17325`, etc.). In E's engine they were derived and therefore
invisible to every solver. They are now free knobs.

## 12. Speed: `fscore.py`
`E2.eqfails` rescans all 39,033 equations per candidate. Only equations touching a nonzero atom
can differ from the constant-only baseline. `fscore` scores in O(#touched eqs): **0.322s -> 0.000s,
2331x**, output verified identical on base and perturbed points. All scans below use it.

## 13. Step 2 — root-firing configurations priced

Deliverable ON-leaves **24601 (A-side / block 178)** and **2081 (B-side / block 21)**: one per
root slot. cfg0's `{1530,1603}` are both B-side (confirmed against `blocks8.json`).

`enum2.py`, from the deliverable base:

    turn OFF 2081                : 38872        turn OFF 24601 : 38909
    turn OFF both                : 38776
    turn ON each of 254 leaves   : best 39000 (leaf 4287), histogram 38961..39000, NONE above baseline
    turn ON best-40 pairs (780)  : best 38975,  NONE above baseline
    simsolve at the base         : 39008

**The deliverable is a strict local maximum: adding and removing both cost.**

## 14. The methodological finding — atom count is the wrong objective
`simsolve` zeroes every bad atom and *loses* 18 points. The deliverable's 8 nonzero atoms produce
only **7 failing equations** because they **cancel inside the equations**. So I solved in
**equation space** (`eqsolve.py`) — choose knob deltas making each equation total zero, allowing
nonzero atoms. Nobody had run this.

    19 knobs, 86 equations; all 5 PIN knobs affine; every equation knob-dependent; 0 unfixable
    exact full solve : FAILS, divisibility on equation 29125  (rhs % -P != 0)
    greedy           : keeps 79/86 -> score 39026

i.e. **the deliverable already is the equation-space optimum over this knob set**, and its 7
failures are structurally forced. Widening to 162 affine knobs / 999 equations (`eqsolve2.py`):
full solve **core infeasible**. The 39,026 barrier is an arithmetic obstruction, not a search
failure — the sharpest characterisation of it so far.

## 15. RETRACTION of my own §9
Last round I wrote *"the root gate never fires at cfg0 or at any configuration E enumerated from
it."* **That is false.** `agentE_work/runs/chanenum.log`: E's 178-channel representatives are
`47` and `112`, both A-side; with `1530`/`1603` B-side already ON, every `mask(1,·,·)` row puts
leaves on both root slots. E enumerated **12 root-firing configurations** (scores 38872-38987).

Corrected reading of E's monotonicity:
- E's own data shows removing leaves also hurts: cfg0 39005 -> cfg5 38917 -> cfg7 38849.
  So the claim is really *"adding channels to a tuned base costs"*, which **does** hold at the
  deliverable base too (best add 39000 < 39026). **Instance property, not a 78-side artifact.**
- The real limitation is the **instrument**: `simsolve` overwrites the tuned handle values, so it
  cannot express a tuned-handle optimum, and tops out near 39005-39008 from *any* base. 39,026 was
  produced by tuned handles. My "78-side artifact" framing was wrong; the coordinator should
  correct FLEET.md.

## 16. Secondary — the slot oracle for F is a NEGATIVE result
`orefine.py`: 136 bases (one saturating each tree stage), corrected engine, signature = delta of
each base's own bad-atom set, ON-corrected pairing rule. Common refinement **32 blocks**
(142,36,16,9,8,8,4,3,2x6,1x18) vs the old 8 — four times finer, **and not valid**:

    partition            crossings vs tree96      vs F's decoded slot pairs
    old 8-block                0                  3 consistent, 0 cut across
    new 32-block              10                  3 consistent, 2 CUT ACROSS

Disagreements: **stage 27596** (F 10|12, oracle 20|2) and **stage 27257** (F 8|3, oracle 9|2).
Stages the oracle splits that are missing from `mux_wiring.json`: `10136`, `15431`, `21279` —
all trivial 1|1 splits of 2-leaf stages, i.e. no new information for the 56 undecoded pairs.

**Conclusion: the residual oracle does not help F's inversion attack.** `blocks8.json` (0
crossings, 3/3 consistent) remains the trustworthy artifact; `oracle_for_F.json` carries both
partitions with validity flags and caveats, `xslot_report.json` the per-stage comparison.

My first crossing test was itself buggy: it flagged a block *containing* a whole stage as a
violation, which is fine for a laminar family. Corrected criterion — a crossing is an overlap that
is neither containment nor reverse containment — gives the numbers above (the buggy test said 78).

## 17. Score attempts this round, all below 39,026
    single-leaf from deliverable base ....... 39000
    leaf pairs (top-40, 780 configs) ........ 38975
    simsolve at deliverable base ............ 39008
    equation-space solve, 19 knobs .......... 39026 (equals baseline, does not exceed)
    equation-space solve, 162 knobs ......... core infeasible
    turning base leaves off ................. 38776-38909

## 18. Step 2 measured from a NEUTRAL base — the curve is unimodal, and the peak is root-firing
`rfenum2.py`, all 2^8 block-subsets, one representative per block, exact re-propagation
(256 configurations in 13 s with `fscore`). Neutral base = the deliverable's handles with both of
its ON leaves turned off (38,776).

    n live blocks : 0      1      2      3      4      5      6      7      8
    best raw      : 38776  38791  38804  38787  38773  38761  38746  38718  38655

The score **rises** from 0 to 2 live blocks and falls after. E's "monotone decreasing in live
channels" is an artifact of measuring only *additions* from a base already sitting at the peak.

    root-firing (A-side + B-side) : best 38804        78-side-only : best 38791
    at exactly 2 live blocks      : root-firing 38804 vs 78-side-only 38774

Best configuration in the entire lattice: `(47, 490)` — 47 in block 178 (A-side), 490 in block 21
(B-side). That is **the deliverable's own slot pattern** (`24601` A-side, `2081` in block 21),
recovered independently from a neutral base. Root-firing wins by 30 points at equal channel count.

Correcting my §15 above: E did enumerate root-firing configurations, but **every one of them also
carried cfg0's two B-side leaves**, so all had >=3 live leaves. **E never priced a clean 2-leaf
root-firing configuration.** Its root-firing rows lost to the 3-live-leaf penalty, not to
root-firing. So the answer is **ARTIFACT** — with a different mechanism than my original §9 gave.

(All raw numbers here are ~38,800 because the handles are untuned at these points; the separate
instrument limitation of `simsolve` in §14 is unchanged.)

## 19. SCOPE — which conclusions are stated against `tree96.json`'s numbering
Agent P's independent parse reports the real tree is **383 law-blocks at depth 9**, not 96 stages
at depth 6; `tree96.json` is a **coarsening**. P independently confirms the 256 leaves, the uniform
law, and the **root split 178 | 78**.

Unaffected, because they never used tree96's node set:
- the 8-block refinement `178|41|21|6|3|3|3|1`; the root split 178|78;
- every score, the gate, `engine2`, the equation-space obstruction, sections 11-14 and 18;
- the set-for-set agreement with `mux_wiring.json` at the root / 19538 / 10649 / 21408 — these
  compare leaf **sets**, which survive any renumbering or further subdivision of stages.

Scope-limited, because they are phrased in tree96's numbering:
- **"Of tree96's 88 non-root stages, exactly 2 are ever split; 86 are never cut"** — true of the
  **coarse** object only. The fine tree contains stages my 339 configurations never separated.
  A limit on scope, not a contradiction, alongside the 8-block resolution limit.
- All stage IDs used above (19538, 10649, 21408, 27596, 27257, 30973, 24533, 10136, 15431, 21279)
  are **tree96 labels**, not fine-tree nodes.
- Therefore §16's negative verdict on the 32-block partition is **downgraded from "spurious" to
  "unvalidated"**: both yardsticks I judged it against (`tree96.json`, `mux_wiring.json`) are
  coarsenings, so its 10 crossings and 2 cut-across verdicts may be the yardstick's error rather
  than the oracle's. `blocks8.json` remains the artifact I would stand behind; the 32-block one
  should be re-judged against the 383-block tree before anyone uses or discards it.

## 20. Still running at hand-off (confirmatory only, both provably capped)
Two tails were still grinding on a machine at load ~27 on 4 cores when I wrote this up:
- `rfenum2.py`'s final block (`simsolve` on the top 6 raw configurations) -> `rfenum2.log`
- `eqsolve2.py`'s greedy row-selection over 999 equations (1500 s budget) -> `eqsolve2.log`

Neither can exceed the baseline: every `simsolve`-based result measured in this round or by E tops
out at 39,005-39,008 because the instrument overwrites the tuned handle values (§14, §18). No
candidate file above 39,026 was produced by anything I ran (`ls M_eq*.json M_rf*.json` -> none).
Whoever picks this up can read the two logs for the confirmatory numbers.

## 21. Round-2 file index
Engine//infrastructure : `engine2.py` (the fix), `fast2.py`, `fscore.py` (2331x scorer),
                         `mcore2.py`, `chan2.py`
Gate and diagnosis     : `gate.py`, `diag1.py`, `diag2.py`, `diag1.pkl`,
                         `val_4287.json`, `val_17378.json` (checker.py-verified)
Step 2                 : `enum2.py` (+`enum2_raw.pkl`, `enum2_pairs.pkl`), `rfenum2.py`
                         (+`rfenum2.pkl`, `rfenum2.log`)
Equation space         : `eqsolve.py`, `eqsolve2.py` (+`eqsolve2.log`)
Oracle for F           : `orefine.py` (+`orefine_sigs.pkl`, `orefine.log`), `ofinal.py`,
                         `xslot.py`, `blocks_corrected.json`, `oracle_for_F.json`,
                         `xslot_report.json`
Superseded             : `rfenum.py` (killed: 65k configs on the slow scorer; `rfenum2.py` replaces it)

---

# LOG_M ROUND 3 — the "obstruction" on equation 29125 does not exist. I withdraw it.

## 22. What I claimed at check-in 8, and why it was wrong
I reported *"the exact full solve fails with a divisibility obstruction on equation 29125
(`rhs % -P != 0`)"* and called it "an arithmetic obstruction, not a search failure".
**That was wrong.** The message was produced by **one elimination ordering inside a badly
overdetermined window** — 86 rows against 19 knobs in `eqsolve.py`, 999 rows against 162 knobs in
`eqsolve2.py`. Neither says anything about the instance. I should have tested the row itself
before naming it an obstruction.

## 23. The window-independent test (`eq29125.py`)
For a single row, integer solvability is exact and needs no window:

    sum_f coef_f * d_f = -s0   solvable over Z   <=>   gcd_f(coef_f)  |  s0

Equation 29125: `issq=False`, 17 atoms, constant part 0, and exactly **one** currently-nonzero
atom — **23617**, coefficient 1 — so `s0 = av[23617]`, a 729-digit integer with `v_P(s0)=0`.

    knob set                     candidates  affine  knobs moving eq29125  gcd  g|s0
    eqsolve 19-knob                  23        19            2              1   YES
    eqsolve2 widened                451       162           12              1   YES
    all cone(eq 29125 atoms)         26        23           12              1   YES

**gcd = 1 in every knob set. Row 29125 is individually solvable.** In all three sets the 5 freed
definer vars are present and affine; `x_28730` is the knob with a direct handle on atom 23617
(`23617 = x_28730 - x_17499*x_9413`).

All seven failing rows tested the same way (`eqsub.py` header):

    eq 12231 gcd 1 | eq 12270 gcd 1 | eq 12350 gcd 1 | eq 14584 gcd 1
    eq 18673 gcd 1 | eq 22044 gcd 40490 | eq 29125 gcd 1        -- all seven divide s0

## 24. Removing the window entirely (`eqsub.py`)
For every subset S of the 7 failures: solve "each row in S = 0", **apply it, re-propagate with
`engine2.forward`, and measure the exact score with `fscore`** — so collateral damage is counted
by measurement, not assumed by a window.

    127 subsets solvable, 0 infeasible, largest solvable subset size 7
    in every case the solver really zeroed its targets (k/k)

    fix eq 12350 -> 10 failures (39023)    fix eq 22044 -> 28 failures (39005)
    fix eq 18673 -> 10 failures (39023)    fix eq 14584 -> 34 failures (38999)
    fix eq 12270 -> 11 failures (39022)    fix eq 29125 -> 34 failures (38999)
    fix eq 12231 -> 18 failures (39015)    fix all 7    -> 44 failures (38989)

    best over all 127 subsets: 39023   < baseline 39026

**There is no obstruction anywhere. Every failing equation is repairable; every repair costs more
than it gains.** 39,026 is a strict local optimum in equation space, and equation 29125 is not
blocked — it is merely tied for the *most expensive* row to repair.

## 25. Why the deliverable's residual is so cheap
The 8 nonzero atoms are visible to exactly the 7 failing equations, nested, all with zero constant:

    12270, 12350, 14584 : all 8          12231 : 23616,23617,36661,36662,36663,36664
    18673 : 36660,36661,36662,36663      22044 : 23616,23617,36664
    29125 : 23617 only

The deliverable places its residual so that only 7 equations ever see it. Equation 29125 sees a
single atom, and the knob that moves that atom is `x_28730` — **one of the 5 vars my fix freed**.
Driving it to zero is exactly what E's orientation did by construction (§11.1); it costs 27 extra
failures. That is the same trade the whole engine defect was hiding.

## 26. Corrected answers
- **Q1 what must divide what:** nothing obstructs; `gcd | s0` holds for all 7 rows.
- **Q2 which knobs move either side:** 12 affine knobs on eq 29125; `x_28730` directly.
- **Q3 eq 29125 or the window:** **the window**.
- **Q4 "core infeasible" at 162 knobs — instance or widening:** **the widening** (999 rows vs 162
  knobs, 6:1 overdetermined). Now proven, not argued: all 127 subsets are feasible.
- **Q5 knob set:** `[642, 7068, 28730, 29854, 31864]` are in it and all affine, checked per run.

## 27. On agent K's claim (unadjudicated, not accepted as input)
K reports the 7 failures are one gate's off-pins forcing the root's two inputs equal so the root
check vanishes identically, and that no configuration can make a stage degenerate, so the freedom
needs broken pins. **My measurements neither confirm nor refute this** — the residual side cannot
see intermediate stage values (§7 of round 1 established that for E's engine and it still holds).
What I can say is only consistency: the residual is confined to 7 equations, and the one repair
that most directly "restores" a check (zeroing atom 23617 via `x_28730`) is among the most
expensive. A real test would need the circuit side: check whether the 8 nonzero atoms are exactly
the off-pins of a single gate, which is F/K's decode to make, not mine.

## 28. Score attempts round 3 — all below or equal to baseline, none above
    every single-row repair ................ 38999 - 39023
    every subset of the 7 (127 of them) .... best 39023
    all 7 simultaneously ................... 38989
No candidate file above 39,026 was produced (`ls M_sub_*.json` -> none).

## 29. Superseded jobs, stopped
`eqsolve2.py`'s greedy row-selection and `rfenum2.py`'s trailing `simsolve` pass were still
grinding when round 3 finished. Both are **superseded and were stopped**:
- `eqsolve2`'s greedy answers the window question badly (it *is* the overdetermined window);
  `eqsub.py` answers it window-free and completely.
- `rfenum2`'s `simsolve` tail is capped at ~39,008 by the instrument limitation of §14/§18.
Their partial logs (`eqsolve2.log`, `rfenum2.log`) are kept; the raw block-lattice table in
`rfenum2.log` (§18) is complete and unaffected.

## 30. Round-3 files
`eq29125.py` (single-row gcd test, 3 knob sets) · `eqsub.py` -> `eqsub.pkl`, `eqsub.log`
(127-subset window-free sweep). Nothing above 39,026 produced. No other agent's directory was
touched; no git commands were run.

---

# LOG_M ROUND 4 — x_7068 identified; alternative placements priced

## 31. x_7068 is a definer that is NOT a handle — collateral, not a fifth corruption
Agent P reports the deliverable corrupts four handle variables `x642, x28730, x29854, x31864`,
exactly the four of 3,707 handles for which `P` does not divide the value. Those are four of my
five freed definer variables. **`x_7068` is the odd one out.** Four independent lines say it is
*collateral damage from `x_642`*, not an independent corruption:

**(a) Definer FORM differs.** The four corrupted handles are defined by a product, or bare:

    x_642    <- x_642   - x_28599 * x_17325       x_29854  <- x_29854 - x_22665 * x_1329
    x_28730  <- x_28730 - x_17499 * x_9413        x_31864  <- x_31864          (bare)
    x_7068   <- x_7068  - x_2099 - 7376877 * x_642        <-- LINEAR COMBINATION

`x_7068` is the only one defined by a linear combination, and that combination **references the
corrupted handle `x_642`**.

**(b) MAGNITUDE is normal.** `x_7068` is **90 digits**, the same scale as `x_2099` (89 digits),
the other term of its own definer. The large corrupted handles are 723-724 digits.

**(c) Its atom is dominated entirely by the x_642 term.**

    atom 23616                                  : 730 digits
    the x_642 term alone, -7376877 * x_642      : 730 digits
    atom 23616 + 7376877*x_642  ( = x_7068 - x_2099 ) :  89 digits

So atom 23616 is nonzero **because `x_642` is corrupted**, not because `x_7068` is.

**(d) The deliverable did not move x_7068 to compensate.** For atom 23616 to be satisfied given
the corrupted `x_642`, `x_7068` would have to be **730 digits**. It is **90**. The deliverable
left `x_7068` at its natural value and let the atom carry the corruption.

**Conclusion: the deliverable corrupts exactly FOUR handles.** My five freed variables are those
four plus one collateral combiner. `P` does not divide `x_7068` either — but that is automatic for
any linear combination containing a corrupted term, and `x_7068` is not product-defined, so it
would not sit in a handle population built from product definers. That P found exactly four
failures among 3,707 handles is **consistent with `x_7068` not being in that population at all**,
which is what its definer form independently says. Two decompositions sharing no data agree on the
same four variables, and the fifth is explained.

(Stated as structure I can verify in my own frame; I cannot read P's handle list and am not
asserting how P's population was built.)

## 32. Alternative placements priced (`engine3.py`, `place.py`) — none beat 39,026
`eqsub` priced *repairs of the current placement*. This prices *different placements* using the
cancellation mechanism the deliverable itself exploits (8 nonzero atoms seen by only 7 equations,
versus E's 2-nonzero-atom state failing 28: more atoms, fewer failures, because they cancel).

**Method.** `engine3.Eng(demote)` generalises engine2 to an arbitrary demotion set. The property
that makes the search sound: **demoting an atom and seeding its variable with its current value
leaves the state bit-identical**, so a demotion is score-neutral and purely adds a degree of
freedom. Verified for every candidate (`0 vars differing, score 39026, 8 bad atoms`); any
candidate failing that check is skipped rather than trusted.

**Candidates.** Atoms that are (i) currently zero, (ii) present in a failing equation, (iii) a
definer (hence demotable): **10 of them** — `11876, 11877, 11878, 11879, 20448, 20451, 20453,
23619, 23623, 36657`. Demoting one frees a new variable that can be moved to *cancel* the existing
nonzero contribution inside that shared equation, instead of driving atoms to zero.

**Result** — for each, equation-target subsets of size 1 and 2 solved, applied, re-propagated and
scored exactly:

    extra=()         39026     extra=(20448,)  39026     extra=(23619,)  39026
    extra=(11876,)   39026     extra=(20451,)  39026     extra=(23623,)  39026
    extra=(11877,)   39026     extra=(20453,)  39026     extra=(36657,)  39026
    extra=(11878,)   39026     extra=(11879,)  39026
    BEST OVERALL 39026  (baseline 39026)   -- "via None": no subset improved on the baseline

**A structural constraint found on the way.** Equations **12270** and **18673** have **zero
demotable zero atoms** — every other atom in them is either already nonzero or not a definer. So
those two failures cannot be addressed by this move at all, whatever values are chosen. Any
placement search that works by adding cancelling atoms inside the failing equations is blocked on
2 of the 7 from the start.

## 33. SCOPE of the placement result — what it does and does not price
State this plainly, because the difference is the whole remaining question.
- `eqsub` (§24) prices **repairs of the current placement**: 127 subsets, all feasible, best 39,023.
- `place.py` prices **alternative placements in a local neighbourhood of the current one** — those
  reachable by demoting ONE additional atom drawn from the current failing equations, with
  equation targets of size <= 2.
- **It does NOT price a genuinely different placement**: corrupting a *different set of handles*
  altogether. That space is not reachable by adding atoms to the current failing equations; it
  requires choosing a different 4-handle corruption, and the handle population is P's object, not
  something I can enumerate from the residual side. **That is the space still untested, and it is
  where anything above 39,026 would have to live.**
- Neither run explored subsets of size >2 for the enlarged knob sets, on the eqsub evidence that
  larger targets are monotonically worse (7 failures -> 44 when all 7 are targeted).

## 34. Round-4 score attempts — all equal to or below baseline
    every single extra demotion x size-1/2 equation targets (11 placements) ... 39026, none above
No candidate file above 39,026 produced (`ls M_place_*.json` -> none).

## 35. Round-4 files
`engine3.py` (configurable demotion set + validate) · `place.py` -> `place.pkl`, `place.log` ·
`placecands.pkl` (the 10 candidates). No other agent's directory touched; no git commands run.

---

# LOG_M ROUND 5 — candidate-agnostic pricer built, calibrated, and used on L's sites

## 36. The primitive (`price.py`, `engine3.py`)
Input is a set of HANDLE variables only. The collateral demotion is **derived**, not supplied:
corrupting handle h frees it, and any variable whose definer atom references h would otherwise
absorb the corruption, so it is demoted too. `closure(handles, depth=1)` computes that.

**Calibration, path 1 (values supplied):** from `[642, 28730, 29854, 31864]` alone,
    closure -> freed [642, 7068, 28730, 29854, 31864]   == engine2's PIN
             -> demote [23616, 23617, 36659, 36663, 36664] == engine2's DEMOTE
    score 39026, fails [12231,12270,12350,14584,18673,22044,29125], 8 bad atoms, 0 vars differing
**The 5th demotion is derived from the 4 handles.** PASSED.

## 37. My first tuner FAILED calibration, and I did not report its numbers
The value-tuning path initially returned **39,008 for the deliverable's own site** — i.e. it could
not rediscover 39,026 from the site that is known to reach it. Twelve candidates all read 39,008.
Per L's caution ("a site priced with unset handles is not a negative result"), those readings were
the tuner's floor and I did not report them as measurements. Diagnosis instead:

`tunediag.py` — and it **refuted my own written hypothesis**. I had assumed the tuner failed
because fixing equations requires accepting newly-broken ones. Measured:

    uncorrupted baseline : 39008, 25 failures, bad atoms [23618,34120,36660,36661,36662]
    deliverable          : 39026,  7 failures
    fixed by deliverable : 18      NEWLY BROKEN by deliverable : 0
    the deliverable's 7 failures are all already failing at baseline

**There is no trade at all** — the deliverable fixes 18 and breaks none. So the target was inside
the searched space and the fault was mine.

`tunediag2.py` — model or solver? Feed the deliverable's own delta (10^724..10^735 on the freed
vars) into the affine model built from +1/+2/+7 probes:

    applying the delta -> 39026, 0 vars differing
    model vs actual on the fixed equations: 12 agree, 0 disagree

**The affine model is exact at a delta of 10^728.** The model was sound; the SOLVER was wrong.
Cause: the knob set had ~40 members, most not globally affine, which let the sparse solver pick
degenerate solutions satisfying the targeted rows while wrecking others.

## 38. Corrected tuner (`tune2.py`) — knobs = the freed handles ONLY
    TUNER CALIBRATION: base(untuned) 39008 -> TUNED 39026   (5 knobs, 18 greedy sols, 4s)  PASSED

The pipeline is now validated end to end: from four handle names, with no values supplied, it
recovers 39,026.

## 39. L's twelve sites, priced with a VALIDATED tuner
    handles                          rows_target   base    TUNED
    [642,28730,29854,31864]  (calib)      25      39008    39026
    [10509,20157,32245,33044]              0      39008    39008
    [9541,19546,25227,31891]               0      39008    39008
    [3260,11588,30400,37248]               0      39008    39008
    [2493,3022,6019,15174]                 0      39008    39008
    [1405,3052,4806,16433]                 0      39008    39008
    [9337,17894,23336,33996]               0      39008    39008
    [19053,21505,22193,23910]              0      39008    39008
    [10074,16399,16800,35694]              0      39008    39008
    [1768,6389,26662,31362]                0      39008    39008
    [6254,7439,21115,38560]                0      39008    39008
    [1079,15006,15333,32131]               0      39008    39008

`rows_target` = how many of the 25 uncorrupted-baseline failing equations the site can even MOVE.

## 40. THE STRUCTURAL RESULT — a hard filter for the candidate generator
**The deliverable's site can move all 25 baseline failing equations. All eleven other sites can
move ZERO of them.** Their 39,008 is therefore not "the tuner found nothing": corrupting at those
sites cannot fix any failing equation, because their atoms do not appear in any of them. The only
possible effect is to ADD failures.

> **A site can help only if its corrupted atoms appear in the equations that fail at the
> uncorrupted baseline.**

The 25 baseline failing equations (relative to the deliverable's free inputs):

    2554, 5324, 6816, 8124, 8680, 9041, 9123, 9421, 11226, 12231, 12270, 12350, 14584,
    15558, 18673, 21000, 22044, 22534, 22997, 28929, 29125, 29330, 32026, 35512, 38051

L's incidence measure is not this quantity: its top 12 by incidence are all 0-incident to these.
Filtering the 378 on this criterion should cut it to a small set, and anything with
`rows_target = 0` can be discarded without pricing.

Caveat: the 25 are relative to THIS baseline (E's orientation from the deliverable's free inputs).
A different free-input configuration would have a different failure set.

## 41. Throughput
    price_given (values supplied) : 0.53 s/candidate  -> ~6,700/hour single-core
    tune        (values tuned)    : 1 s for a 0-incidence site, 4 s for a fully incident one
                                    -> ~900-2,700/hour single-core, 4 cores available
So L can emit **thousands**; list size is not the constraint. Better: L can pre-filter on the 25
equations above and emit only incident sites.

## 42. Round-5 files and score
`price.py` (Pricer/TunedPricer, closure) · `pricetest.py` · `pricerun.py` (+`pricerun.pkl`) ·
`tunediag.py`, `tunediag2.py` · `tune2.py` (+`tune2.pkl`, `tune2.log`).
Nothing above 39,026 produced; no `M_site_*.json` written. Baseline stands.

---

# LOG_M ROUND 6 — the filter is baseline-independent, and the incident set is 0.28%

## 43. THE CAVEAT RESOLVED: the two baselines are IDENTICAL (`basecmp.py`)
The 25 were computed in E's orientation, which is the known-defective frame. Recomputed
independently against the deliverable's own baseline:

- **Baseline A** — E's ORIGINAL orientation, full forward from the deliverable's free inputs.
  Every definer atom forced to zero, including the 5 the deliverable needs nonzero.
- **Baseline B** — start from the deliverable's actual vector in the CORRECTED engine and
  **un-corrupt it in place**: set each freed variable back to the value its own definer atom
  prescribes, iterated to a fixpoint (needed because x_7068's definer references x_642), then
  re-propagate. Cofactors and all else keep the deliverable's values.

```
BASELINE A: score 39008, 25 failures, 5 bad atoms
BASELINE B: score 39008, 25 failures, 5 bad atoms
|A n B| = 25    in A not B: []    in B not A: []    IDENTICAL: True
```
(Baseline B also confirms the un-corruption worked: all 5 demoted atoms return to zero, and
x_31864 returns to 0, correct for a bare definer.)

**The filter is not a property of E's orientation.** L can use it. `baseline_sets.json`.

**Remaining caveat, unchanged and still real:** both baselines share the deliverable's FREE
INPUTS. The 25 are a property of the instance *at this free-input configuration*, now verified
across two orientations. A materially different free-input configuration could give a different
set; that is untested.

## 44. How thin is the incident set? 32 of 11,307 (`incid.py`)
Characterised from the equation side only -- no site enumeration, so the boundary with L holds.
Moving a freed handle u changes every atom that mentions u, so

    u is incident  <=>  occ[u]  intersects  {atoms appearing in the 25 baseline failures}

    atoms appearing in the 25 equations ......... 78
    variables touching at least one such atom ... 131
    product-defined variables (handle shape) .... 10381   bare-defined 926   total 11307
    *** INCIDENT HANDLE POOL ................... 32  (0.28%) ***
    the deliverable's four 642,28730,29854,31864 : ALL FOUR IN THE POOL

    [642, 1627, 1844, 1956, 2218, 2892, 4863, 6480, 7062, 7945, 9629, 10861, 11425, 15422,
     16495, 21279, 21718, 23538, 23642, 23754, 23822, 26732, 28098, 28730, 29305, 29854,
     30175, 31465, 31864, 33001, 35619, 37413]

### Cross-check against the measured data — passes
    L's 11 sites, 44 distinct handles : 0 of 44 in the pool  <-> measured rows_target = 0 for all 11
    the deliverable's 4               : 4 of 4 in the pool   <-> measured rows_target = 25
The criterion derived from the equation side agrees with the pricing measurement on every site
tried, from both directions.

## 45. What this does and does not say
- **Does:** the search space collapses from C(11307,4) = 6.81e14 four-handle sites to
  C(32,4) = **35,960** -- a reduction of **1.9e10**. At the measured 1-4 s per tuned site that is
  **~20 core-hours, ~5 hours on 4 cores: the alternative-placement space is now EXHAUSTIVELY
  PRICEABLE.** That is a far sharper statement than "no placement below 7 found".
- **Does not:** incidence does NOT force the deliverable's site. It is one of 35,960, not unique.
  What is established is that 99.72% of handles cannot participate at all.
- **Necessary, not sufficient:** `occ[u] n A25 != {}` is a necessary condition, so 32 is an
  UPPER BOUND on the helpful pool; the sufficient test is the measured `rows_target`. Some of
  the 32 will price out at 0 rows once tuned.

## 46. Round-6 files
`basecmp.py` -> `baseline_sets.json`, `basecmp.log` · `incid.py` -> `incident_pool.json`,
`incid.log`. No candidate generator built. Nothing above 39,026. No other agent's directory
touched; no git commands run.

---

# LOG_M ROUND 7 — baseline gap explained; my 32-pool corrected to 103; the lead priced

## 47. THE BASELINE DISCREPANCY — explained, and it is ONE atom (`lcrit.py`)
My baseline fails 25; L's fails 13, a strict subset. Constructing L's baseline in my own frame
(the deliverable with its 16 tuned handle/cofactor variables zeroed):

    L-style baseline (16 zeroed) : score 39020, 13 failures, bad atoms [23616, 23618, 36660, 36662]
    my baseline (un-corrupt)     : score 39008, 25 failures, bad atoms [23618, 34120, 36660, 36661, 36662]
    L-style is a strict subset of mine; in mine not L's = 12, in L's not mine = 0

**Every one of the extra 12 fails for exactly one reason: atom 34120.** Checked per equation --
in my baseline each of the 12 has 34120 as its only nonzero atom; in L's, none of them has any.

    eq 5324, 9041, 11226, 15558, 21000, 22534, 22997, 28929, 29330, 32026, 35512, 38051
      -> nonzero atoms in MY baseline: [34120]   in L-style: []

Atom 34120 is one of the two atoms `x_7068` touches (the other is 23616; measured at
check-in 30). **My un-corruption puts `x_7068` back on its definition
`x_2099 + 7376877*x_642` -- a 735-digit value -- which makes 34120 nonzero. L's zeroing leaves
`x_7068` at the deliverable's own 90-digit value, where 34120 is zero.**

> **Answer to the question asked: un-corruption propagates further than zeroing.** The extra 12
> are NOT an artifact of my method and L's zeroing is not leaving something wrongly satisfied.
> In the honest uncorrupted machine `x_642` is large and `x_7068` follows its definition, so
> atom 34120 genuinely is nonzero and those 12 genuinely do fail. **L's decision to filter on
> the union (25) was right, and is now justified rather than merely cautious**: a site that can
> only fix one of those 12 would be wrongly discarded by the 13-set.

### Corollary — the deliverable's 18 fixes decompose exactly
    the 18 it fixes = the 12 above (all killed by zeroing atom 34120 via x_7068)
                    + 6 more [2554, 6816, 8124, 8680, 9123, 9421] (via the handle corruptions)
That is why the deliverable holds `x_7068` at a small value instead of its definition: **12 of
its 18 fixes are bought with that one variable.** It also explains why `x_7068` was the single
collateral demotion -- it is not incidental, it is doing most of the work.

## 48. CORRECTION TO MY OWN 32-POOL: the real pool is 103, not 32
My "32 incident handles of 11,307 (0.28%)" restricted the population to **product/bare-defined**
variables -- the identical blind spot I diagnosed for `x_7068` at check-in 18, now at scale.
Recomputed over **all definer forms**, and counting both routes by which a handle reaches a
target equation (its definer atom being in one, or an atom containing it being in one):

    CORRECTED POOL: 103 incident handles of 30,383 definer variables (0.34%)
    all 32 of my old pool survive; 71 handles were MISSING from it
    deliverable's four: all present     stage checks 23754/35619/9629: all present

Exact-by-definer-atom counts, for the record: 59 incident handles against the 25-equation
baseline, 22 against L's 13-equation baseline (the 13-set is a subset of the 25-set). All six
handles the coordinator relayed from L appear in my 25-set derivation, independently obtained.

**The consequence for the enumeration I proposed: C(32,4) was the wrong space.** It was a
shape-restricted subset, not an upper bound as I claimed. I withdraw the "0.28%" figure.

## 49. Fast tuner (`sweep.py`) -- validated, ~10x faster
Scoring a greedy prefix used a full forward + badatoms; only the freed variables move, so an
incremental `resid_delta` gives the same bad-atom dict. The definer-level user map is built once
globally instead of per site. Validated: calibration `39008 -> 39026` in **0.4 s**, and a full
re-propagation of the same seed returns **39026**, so incremental scoring is exact.

## 50. THE LEAD, PRICED: all 98 five-handle supersets -- nothing above 39,026
Deliverable's four plus one handle from the pool, stage checks first in the ordered prefix:

    +x23754  rt 20  rows 25  -> 39026     <-- stage check
    +x35619  rt 17  rows 25  -> 39026     <-- stage check
    +x9629   rt 16  rows 25  -> 39026     <-- stage check

    DISTRIBUTION over 98 priced five-handle sites
      39026: 89      39012: 1      39011: 8
      above 39026: 0        priced out at 0 rows once tuned: 0 of 98

**The three stage checks come back exactly at baseline.** That is consistent with L's correction
that they are vacuous at the deliverable's configuration (`sel_ab(x27994)=0`): a vacuous atom
supplies a free additive term, and the tuner already reaches 39,026 without it, so it adds
nothing. Adding freedom is score-neutral at best (89 cases) and sometimes worse (9 cases,
39,011-39,012) where the extra demotion perturbs the greedy path. **Zero sites priced out at 0
rows** -- every handle in the pool is genuinely incident, so the pool is not padded.

## 51. O's lattice target PRICED — the region can be made to hold, collateral is 44-49
I do not have O's integer shift vector, so I re-derived the target in my own frame: solve the
region equations exactly with a knob set that **includes the free carriers**, which is precisely
what my earlier `eqsub` run lacked (it used only the freed handles and reported 38,989).

Confirmations of O's setup, independently in my frame:

    a23616 = x_7068 - x_2099 - 7376877*x_642      a36660 = 5113045*(x_7075*x_9118) - x_29854
    a23618 = x_4432 - x_19964 - x_28730           a36662 = x_7075 * x_8731
    carriers x_8731, x_9118, x_4432 are all FREE and all AFFINE here
      x_8731 -> atom 36662 only     x_9118 -> atom 36660 only     x_4432 -> atoms 8721, 23618
    the deliverable's 7 failures are a strict SUBSET of L's 13-equation region
    x_17499 == p EXACTLY  (corroborates O's account of the eq29125 elimination)

**Result — the target is reachable, and it is not free:**

    target: all 7 currently-failing hold  -> SOLVED 7/7   score 38989  (44 equations fail)
    target: all 13 region equations hold  -> SOLVED 13/13 score 38984  (49 equations fail)
    every subset of the 7, carriers included -> best 39026 (nothing above baseline)

So **O's inversion is confirmed from the residual side: the whole region CAN be made to hold
simultaneously.** The collateral, which O's model cannot express, is **44 equations** for the
7-target and **49** for the full 13-target. Net against the 39,026 baseline: **-37 and -42.**

### The caveat that matters, and it is not a quibble
**I priced *a* solution to O's target, not necessarily O's δ₀.** The solution set is a lattice
coset and my solver picked a different point: my shifts are **~4,200-4,558 bits**, O's are
**2,419-2,440 bits** — roughly half the size. Collateral plausibly grows with shift magnitude
(a larger shift perturbs more downstream atoms), so **O's smaller δ₀ could price materially
better than mine, and it cannot be inferred from bit-sizes alone.** To price O's actual point I
need its integer vector. That is the single highest-value thing to relay to me.

Knobs my solution moved (11 for the 7-target, 17 for the 13-target) include the cofactors
x_950, x_1329, x_3629, x_6418, x_6947, x_8976 and the carrier x_9118 — i.e. the solve does reach
for the zero-collateral carriers, but not exclusively.

## 52. Round-7 files
`lcrit.py` -> `lcrit.json` (baseline reconciliation) · `exact15.py` -> `exact_incident.json` ·
`sweep.py` (fast validated tuner) · `pricelead.py` -> `pricelead.pkl` (98 five-handle sites) ·
`pricelead2.py` (six-handle, stopped for the O target) · `pricedelta.py` -> `pricedelta.log`.
Nothing above 39,026 produced anywhere. Baseline stands.

---

# LOG_M ROUND 8 — O's delta0 priced. Minimal representatives do NOT reduce collateral.

## 53. What I confirmed of O's model (empirically, in my frame)
The right instrument is a +1 probe and re-propagation, not the `occ[]` listing -- `occ[]`
includes atoms whose coefficient is zero at this configuration (e.g. `x_21279 = 0` kills the
apparent `x_9118 -> x_25297 -> x_2099 -> a23616` path). Measured:

    move x_9118 by +1 -> changes region atom 36660 ONLY.  outside atoms changed: NONE
    move x_8731 by +1 -> changes region atom 36662 ONLY.  outside atoms changed: NONE
    move x_7068 by +1 -> changes region atom 23616  AND outside atom 34120
    move x_4432 by +1 -> changes region atom 23618  AND outside atom 8721

**O's zero-collateral claim for the two free carriers is CONFIRMED.** The two derived carriers
leak to exactly **one** outside atom each. (I first read `occ[]` and wrote that none of the four
was private; that was wrong and the probe corrects it.)

**And atom 34120 is the one that drives 12 of the 25 baseline equations (LOG 47).** So `x_7068`
is doing double duty: it buys 12 of the deliverable's 18 fixes, and it is also the only route
for the a23616 shift. That is the tension in one variable.

## 54. delta0 priced -- it does not land, and the reduction does not help
Applied from the deliverable in my frame, carriers shifted and O's `z` on the private set:

    variant                                   score   region still failing
    A  direct, z as increments                38998   7/13
    B  direct, z absolute                     38992   13/13
    C  MINIMAL REPRESENTATIVES, increments    38993   12/13
    C  minimal representatives, absolute      38992   13/13
    D  reduce only the x_7068 carrier         38995   10/13
    E  reduced carriers only, no z            38992   13/13
    E  free carriers only                     39022   11/13

Reduction achieved: **d1 2440 -> 22 bits** (`-3228258`, mod 7376877), d2 2419 -> 255,
d3 2406 -> 255, d4 2429 -> 253, with the compensations that hold the region atoms fixed
(`x_642 += (d1r-d1)/7376877`, `x_28730 -= (d2-d2r)`, `x_9413 += (d2-d2r)/p`).

> **The bet is refuted. Shrinking the first carrier from 2440 bits to 22 did not reduce
> collateral -- it made it slightly worse (38,993 vs 38,998). Collateral here is structural,
> not magnitude-driven.** The cost is that `x_7068` moves atom 34120 at all, and 34120 costs
> 12 equations whether it is perturbed by 2^2440 or by 2^22.

I also swept all 12 sign/mode conventions (carrier +/-, z +/-, z as increment/absolute/omitted):
**best remains 38,998 with 7 of the 13 region equations still failing.** Under no interpretation
I tried does delta0 make the region hold in the realisable system.

## 55. But the target itself IS reachable -- I reached it independently
My own region solve (LOG 51, `pricedelta.py`, carriers in the knob set) makes **13/13 hold**, at
score **38,984** (49 equations of collateral). So the discrepancy is not "the region cannot be
satisfied"; it is that **O's particular delta0 does not satisfy it once propagation is real**,
while a different point in the same coset does, expensively.

Two readings, and I cannot separate them from my side alone:
  (a) my mapping of O's boundary shifts onto variable increments is not the one O intends --
      O should state whether the a23616 shift moves `x_7068` or `x_2099`, and whether `z` is an
      increment from the witness or an absolute value; or
  (b) O's 9-atom/13-equation region model omits a constraint that full propagation enforces.
The probe results in section 53 are evidence for (b) being small -- the carriers behave exactly
as O describes -- which makes (a) the more likely gap.

## 56. Standing summary of measured prices (nothing has beaten 39,026)
    deliverable                                  39026
    O delta0, best of 12 interpretations         38998
    my own solve making all 13 region eqs hold   38984
    all 7 failures forced to hold (my solve)     38989
    98 five-handle supersets, best               39026 (89 of 98 exactly equal)

## 57. Round-8 files
`delta0.py` -> `delta0.log`. Nothing above 39,026 was produced, so no assignment was written to
disk and there was nothing to verify with `checker.py`.
