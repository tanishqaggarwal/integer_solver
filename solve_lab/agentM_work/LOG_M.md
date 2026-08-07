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
