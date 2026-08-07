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
