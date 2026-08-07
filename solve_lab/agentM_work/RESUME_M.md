# RESUME_M — agent M. Angle: is E's channel partition the same object as F's tree?

## 0. Verdict (short)
**CONFIRMED.** The residual-side channel partition and the circuit-side subtree partition are
two views of one object. E's untested prediction — that the 16 inert booleans become live at a
configuration that fires their stage — is now **measured true**, and the mechanism is identified.

## 1. Baseline re-verified by me
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
Atom-level: exactly **8 nonzero atoms** `{23616,23617,36659,36660,36661,36662,36663,36664}`.
Its ON-leaf set is `{24601, 2081}`. **I did not beat 39,026.**

## 2. What I ran (all in `agentM_work/`, absolute paths, nothing else touched)
- `mcore.py`   — channel measurement at ANY seed (E's `channels()` generalized; full 5-coord
                 signature over ROWS `[7389,10187,20212,20215,28647]`, not just E's 2-coord projection).
- `xcompare.py`— decompose a class into maximal `tree96.json` nodes; count tree nodes split across classes.
- `sweep1.py`  — 42 configurations, partition + crossing count -> `sweep1.pkl`.
- `refine.py`  — 257 configurations (all-off + each of the 256 leaves alone ON); common refinement
                 of all observed partitions -> `refine256.pkl`.
- `refine2.py` — 297 configurations, CORRECTED refinement (two leaves together iff their labels
                 agree in every config where NEITHER is ON) -> `blocks8.json`. The naive version
                 (`refine.py`) is contaminated: each leaf is uniquely ON in one config, so it
                 splits all 256. Do not reuse `refine256.pkl`.
- `fullsig.py` — wider signature (full bad-atom delta) — too fine, splits all 256 (private pins).
- `ancestors.py` — stage-wire signature — FAILS in E's engine (see 7).
- `mchan.py`   — E's `channels.py` with absolute paths (read-only copy; E's dir untouched).
- `pairscan.py`, `scan_single.py` — score attempts, both below baseline (see 7).
`orient.pkl` regenerated in my dir is **byte-identical** to E's, so all measurements are in E's model.

## 3. The correction that unlocked it
E reported "16 booleans moving nothing at all" at cfg0. Reproduced exactly (178/41/21 channels,
partition identical to `chan_cfg0.json`), but the 16 are **14 inert + 2 already ON**
(`x_1530`, `x_1603`) — E's measurement skips ON bits (`if v0[f]==1: continue`), so they were
counted as "not in a channel" rather than "already firing".

`x_1603` sits in tree node **19538** (9 leaves); the other 8 leaves of 19538 are the inert ones.
`x_1530` sits in tree node **10649** (4 leaves); the other 3 are inert.
That is the exact signature of E's own saturation law: **one leaf ON saturates its whole subtree,
so every other leaf of that subtree measures as inert.**

## 4. The test E asked for — RESULT: they become live
Set `x_1530 = x_1603 = 0` in E's `triple8_seed.json` and re-measure:
| config | 16-group splits into |
|---|---|
| cfg0 (1530,1603 ON) | 14 INERT + 2 ON |
| 1530=1603=0 | **ch9 = node 19538 exactly, ch4 = node 10649 exactly, +3 into a ch24** |
| deliverable config `{24601,2081}` ON | **one single live channel of exactly 16** |

At the deliverable's configuration the 16 are `19538(9) ∪ 10649(4) ∪ 19327(1) ∪ {17760,21074}`
= one live channel. So `178 | 41 + 21 + 16` is a real partition of the 78-side, not numerology.

## 5. The partition IS the tree — and F's circuit confirms it set-for-set
Slot leaf supports computed from **F's `mux_wiring.json` alone** (selector vars -> defining atom's
cone -> intersect the 256 leaves), never from my measurement:

| stage | circuit inA/inB | residual-side block |
|---|---|---|
| **15298** root | 178 / 78 | inA **== block 178 set-for-set**; inB == union of my other 7 blocks |
| **19538** | 6 / 3 | **== blocks 6 and 3, set-for-set** |
| **10649** | 3 / 1 | **== blocks 3 and 1, set-for-set** |
| **21408** | 10 / 11 | at deliverable cfg: live class of 11 **== inB exactly**; inert-inside-21408 **== inA \ {2081} exactly** |

Maximal common refinement over 297 configurations (`blocks8.json`):
**178 | 41 | 21 | 6 | 3 | 3 | 3 | 1**, and 41+21+6+3+3+3+1 = **78**. Of tree96's 88 non-root
stages, exactly **2** are split by it (19538 -> 6|3, 10649 -> 3|1) and both splits are the stage's
own two input slots. **86 stages are never cut; zero arbitrary crossings.**

## 5b. Raw crossing statistics
Intersecting E's cfg0 classes with all 96 stages of `tree96.json`: **every one of the 96 stages
lies entirely inside one class.** Zero stages are split (after removing ON leaves, which cannot
be flipped and are a measurement artifact, not a crossing).

Over the 42 configurations of `sweep1.pkl`, the ONLY tree nodes ever split are **19538 (always 6|3)**
and **10649 (always 3|1)** — and those splits are the nodes' own two input slots, i.e. a
*refinement below* `tree96`'s granularity, not a violation. Same at the deliverable's config,
where node **21408 (21)** splits **11|10** with the ON leaf `2081` inside the 10 side.
**Zero arbitrary crossings in 42+257 configurations.**

## 6. Mechanism (the statement to quote, with knob set and configuration)
Knob set: the **256 boolean leaves** `mcore.bools()` (the booleans of the cone of ROWS
`[7389,10187,20212,20215,28647]`), each probed 0->1 singly.
Statement: at any selector configuration C, two leaves are in the same channel iff their lowest
common **firing** ancestor is the same; equivalently

> **the channel partition at C = the leaf-support partition induced by the tree's slot structure,
> cut at the deepest gate that C saturates.**

Turning one leaf ON inside a subtree S saturates S and splits S's class into S's two slot supports.
This makes the channel measurement (0.6 s for all 256 leaves at once) a **direct oracle for the
binary slot structure of the tree** — including the 56 slot pairs F left undecoded.

## 7. What this does NOT show
- It does not give a score above 39,026. E's `simsolve` run from the deliverable's configuration
  returns **39,008** (its knob discovery is tuned to cfg0's cluster; the deliverable's residual
  atoms `{23618,34120,36660,36661,36662}` under E's orientation are a different cluster).
- E's `forward` from the deliverable's free inputs gives **39,008, not 39,026**: E's orientation
  zeroes the 8 atoms the deliverable leaves nonzero and pushes the defect elsewhere. 23 derived
  vars differ. So E's engine is NOT a superset of the deliverable's configuration. This is the
  most important negative result here.
- Score attempts, all below baseline: `scan_single.py` (each leaf alone ON through `simsolve`)
  38,842; `pairscan.py` (deliverable's 30 handle settings kept, leaf pair swapped over 14x78)
  best 38,944 at `(21266,2081)` vs 39,008 for the deliverable's own pair — handles are tuned to
  the leaf pair, so swapping leaves without re-solving handles is a dead instrument.
- **Resolution limit of the residual oracle: 8 blocks.** It cannot see inside the 178, inside the
  41, or inside 21408 beyond one slot split. Widening the signature to the full bad-atom delta
  support splits all 256 (private leaf pin atoms) — `fullsig.py`.
- **The stage-wire oracle does not work in E's engine.** Recording which of the 96 stages' wires
  change: at all-off 219/256 leaves change NO stage wire; at cfg0 229 change exactly one (the
  root). E's forward is a propagation-with-defect that does not realize intermediate stage values.
  So the residual side alone can never decode a-side internal wiring; F's decode is not replaceable.

## 8. Highest-value next experiment
E's entire channel enumeration (`chanenum.py`, "the empty set wins at 39,005", monotone) was
anchored at cfg0, **whose ON-set `{1530,1603}` is entirely on the 78-side** — the root gate never
fires there. The deliverable's ON-set has one leaf per root slot and scores 21 higher. E's
monotonicity result therefore does not price the configurations that matter.
Next: re-run the channel/representative enumeration from a base with the root gate FIRING, and
fix E's orientation so it can represent the deliverable's 8-atom residual (see §7).

## 9. Files
`mcore.py xcompare.py sweep1.py refine.py refine2.py fullsig.py ancestors.py mchan.py
pairscan.py scan_single.py cfg0_M.json deliv_seed.json blocks8.json blocks2.pkl sweep1.pkl
refine2.pkl anc_alloff.pkl pairscan.pkl orient.pkl LOG_M.md`
Full narrative with every measurement: **`LOG_M.md`**.
