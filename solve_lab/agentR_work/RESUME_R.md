# RESUME_R — agent R. Self-contained. Read §0, then §A; §B is withdrawn, §C is how I got burned.

## 0. Score and standing
- Baseline re-verified by me from cold, twice:
  `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
- **I did not beat it.** Nothing I produced scores above 39,013 on `checker.py`.
- **No infeasibility is claimed anywhere in this directory.** A solution exists.
- My angle produced **two withdrawn results and one durable explanation**. The sections below are
  separated so the next reader can tell which is which without reading the narrative.
- Full narrative with every number and every dead end: `LOG.md`. Raw data: `runs/`.

---
# §A. SURVIVES — pure equation-incidence, no group model anywhere in it

These are facts about `EQUATIONS.txt` as a set of equations over atoms. They do not depend on the
accumulator model that died in §B, and nothing in §B is needed to state or check them.

## A1. The footprint-cost framework
For a set S of atoms held nonzero:
`cost(S) = |equations touched by S| − max #{touched equations a nonzero value vector annihilates}`.
`cost(S)` is the **failing floor** for any assignment whose nonzero-atom set is exactly S.
Beating 39,026 requires a support with `cost(S) ≤ 6`.
Agent F's evaluator model: **9,032 atoms, 39,033 equations, 121,261 incidences.**

## A2. THE DURABLE RESULT — why every configuration-first search bottoms out at 7
Single-atom footprint cost distribution (cost -> number of atoms):
2->1, 3->1, 4->1, 5->1, 6->2, 7->8, 8->36, 9->113, 10->354, 11->863, 12->1481, 13->1793,
14->1835, 15->1327, 16->721, 17->339, 18->112, 19->33, 20->9, 21->1.

Split by kind, which is the whole point:

| | minimum single-atom cost |
|---|---|
| **relational** atoms (≥2 variables — the only ones that can carry the defect) | **7** |
| **boolean-ness** atoms `x·(1−x)` (one variable; assert `x ∈ {0,1}`) | 2 |

> **Every atom with footprint cost ≤ 6 is a boolean-ness atom on a single variable.**
> The cheapest *relational* atoms cost exactly 7 — **3131, 3138, 3588, 8749, 8777** — and
> **atom 3131 is one of the deliverable's own four live atoms.**

**So the deliverable sits at the single-atom relational optimum, and 7 failing is that optimum.**
`|S| = 2` ties but never beats: over the 60 cheapest relational atoms, all 34 equation-sharing
pairs floor at ≥ 7, the best (3131+3138, 3138+8749, 3588+5558, 8777+2569) tying at exactly 7
(9 touched, 2 killed) — cancellation buys back exactly what the larger support costs.

This is the best explanation anyone has produced for why five independent configuration-first
searches all stopped at 7. Reproduce with `footprints.py` then `rank.py`.

## A3. The disconnected-cheapness trap
The cheapest atom of all is **8508 = `x29570·(1−x29570)`, cost 2** — and `x29570` occurs in **no
other atom**. It is cheap *because it is disconnected*, and a disconnected variable cannot carry
the defect (the defect links the leaf selection to the target and must flow through atoms that
connect them). Same for 7888 (`x27026`, cost 5), 8512, 8764 — all zero-collateral.
**Cheapness and load-bearingness anti-correlate here.** This eliminates the four cheapest entries
and is the ranking's own version of "incidence does not price cost". `reach2.py`.

## A4. The 39,029 floor — the arithmetic, not the reachability
Exhaustive over all 253 selector boolean-ness atoms that have a ladder position:
`{x24267 (atom 7887), x33095 (atom 8509)}` have a **4-equation union -> floor 39,029**; three
further pairs floor at 39,027; four tie at 39,026. **These are floors, i.e. lower bounds on
failing equations for those supports. The argument that any of them is REACHABLE is withdrawn
(§B2).** `obstruct.py`, `runs/obstruct.json`.

## A5. Cancellation pricing, stated at measured strength
`gs2.solve`'s single-bit footprint is 3 live atoms in **20 equations with 0 cancelling**; the
deliverable's is **13 equations with 6 cancelling**. Of the 20: **0** can never cancel (all have a
dead partner atom available), but **0** have a partner touching ≤3 other equations. The cheapest
partner anywhere is atom 7954, present in **10** other equations; then 4490/4497 at 11, 4496/4561
at 12, 8259 at 13, 4500/4331 at 14.

**Stated as measured: buying one cancellation TOUCHES ≥10 further equations. Touched is not
failed** — they may cancel downstream, and this lab has repeatedly shown incidence does not price
cost. This is a price floor in *touches* with a ~10x margin against the 13-equation deficit, **not**
a proof the purchase cannot pay. **The margin is what makes the conclusion likely, not the
arithmetic.** `cancel.py`.

## A6. TOOLING FINDING — agent F's evaluator `E` mis-scores the deliverable
| evaluator | failing on `best/new_instance_partial_39026.json` |
|---|---|
| `checker.py` (ground truth) | **7** — 12231, 12270, 12350, 14584, 18673, 22044, 29125 |
| agent F's `E` | **13** — those 7 **plus 2554, 6816, 8124, 8680, 9123, 9421** |

`E` scores the deliverable **39,020 instead of 39,026**. The over-report is **assignment-dependent,
not a constant offset**: re-scoring with `checker.py` (`rescore.py`), `E` and `checker.py` agree
**exactly** on all five `gs2` assignments I generated (39,013 / 39,013 / 39,013 / 39,005, over by 0).
So scores computed with `E` are usable, but **the deliverable's own footprint cannot be read off
`E`**, and anything compared against it via `E` is off by 6. Agent L found the same 13 independently.

## A7. The pins -> ladder lookup, corrected and validated by wire contents
`agentF_work/pins.json` stores each pin as `[[wire, val], [wire, val]]` where the **second** entry
is the x-coordinate, **unreduced**, in the unshifted frame. So
`ladder index = LX[(val_second + S) mod P]`. Validated: **x24601 -> 72, x2081 -> 235**, matching
four agents' independent readings.

Stronger, and measured against the deliverable's actual wire contents (`validate_A2/A3.py`):
* exactly **2 of 256** pin variables are nonzero in the deliverable — x2081 and x24601;
* the **4 coordinate wires those two pins name are all set, all holding exactly the named value**;
  the other 505 coordinate wires are unset;
* exactly **2 of 256 ladder points** have a coordinate on a wire — leaves **72 and 235**;
* the target's coordinates are on wires (x on 13682/22162/24468, y on 10156/18956/30213).

This is a far stronger check of the ON-set than my earlier `fold(k) != T` test, which was
**nearly vacuous** (almost any wrong model also yields `fold != T`). Note this is consistent with
T's and Q's finding rather than in conflict with it: the coordinate wires are free variables the
deliverable assigns, so forcing the selectors to 0 does not clear them, and a leaf pin
`sel·(w − C) − z` lands its coordinate only once `z` is separately forced.

## A8. Solver measurements — survive as measurements of TOOLS, with one dependency flagged
Full tables in `LOG.md` §4. Headline numbers: z3 QF_BV needs **119 s** on an 8-bit-prime, 7-stage
instance with **every selector pinned to the known solution**, and times out with them free;
CaDiCaL 1.9.5 needs **165,725 conflicts** to verify the same pinned instance; CP-SAT (much the best
of the three) solves 8- and 10-bit free instances then fails at 12, 14, 16, 20, 24, 28, and returns
MODEL_INVALID at 31 bits because it cannot represent the arithmetic; bit-blasted CNF grows as
`≈863·m²` clauses per stage.

**Dependency to flag:** these ran on `sibling.py` instances built to my chain model, which §B
refutes as a model of the circuit. They therefore measure *solvers on modular chord-law chains*,
which is an honest measurement of solver capability but not automatically a measurement of this
instance. The clause-count extrapolation (≈1.4x10^10 clauses, ≈600 GB of DIMACS against 29 GB of
disk) is the part least exposed — bit-blasting 256-bit modular multiplication is quadratic
whatever the topology — but making it model-free would mean recounting multiplications directly
from `EQUATIONS.txt`. I did not do that. Treat A8 as suggestive, not as A2-grade.

---
# §B. WITHDRAWN — died with the accumulator model. Do not carry forward.

## B1. What died and why
`solve2.py` seeded the accumulator chain at `L_0`. **Leaf 0 is absent from the deliverable while
both live leaves are present** — if the fold were seeded at `L_0` the deliverable's live-leaf set
would be {0, 72, 235} and the ON-set would carry a `2^0` term. It does not; four agents read it as
`2^72 + 2^235`. **`A = L_0` is refuted.**

Worse for the model generally: **no accumulator value of any kind appears on any wire.** Not the
fold `L72 + L235`, not `L0 + L72`, not `L0 + L72 + L235`, in any frame. **The deliverable holds the
inputs and the target and nothing between them** — consistent with T's and Q's finding that
routing is a constraint that is not propagated, and with K's finding that root slots are not pinned.
Three agents' models failed this way in one round; mine is one of them.

## B2. The specific claims withdrawn
* **The four `(t1, t2)` root results** of `LOG.md` §16.5, including the 39,029 pair. Roots of a
  system built on the refuted seed. **Nothing was routed to agent M and nothing should be.**
* **The degree-collapse argument** (`collapse.py`) *as an instance claim*. Still true of the
  siblings, but the siblings are my own construction and the collapse depends on the same model.
* **"Relaxing a selector leaves the mux atoms satisfiable."** Rests on the mux form
  `acc' = acc + b·(S − acc)`, which is model, not measurement.
* **The 2^gap backward-solve obstruction** (`tradeoff.py`). I had already falsified this one myself
  before the model died — it assumed the intervening selectors are arbitrary when I choose them.
  Retained in `LOG.md` §16.3 only because the reasoning error is the instructive part.

## B3. Point-level identities — true as identities, NOT as circuit semantics
These I measured on the point set, and they remain true of those points. What does **not** follow,
and what I withdraw, is the inference that *the circuit computes with them*:
* `X = x + K/3` removes the offset; the stage law becomes the plain chord law (200 random pairs).
* All 248 forced pin points and the target satisfy `Y² = X³ + B` with
  `B = 64019533680030876408443198762210829058751700634554282185987325820393598524794`
  (fitted from 2, verified on 246 + target, 0 exceptions).
* That law is commutative and associative (200 random triples each).
* The 248 points form one doubling chain `L_i = 2^i·L_0`, recovered as 9 doubling-closed pieces
  (111,61,28,12,11,9,8,5,3) with 8 two-doubling splices -> `ladder.json`, 256 points.
* Cornacchia on `4p = L² + 27M²` gives order
  `N = 115792089237316195423570985008687907852837564279074904382605163141518161494337`,
  256-bit prime.
* Exhaustive Hamming weight ≤ 6 meet-in-the-middle: **no `k` of weight ≤ 6 satisfies `k·L_0 = T`
  in that group** (108 s). This stands as a statement about the group, **not** about the circuit —
  same withdrawal class as agent Q's six search programs.

---
# §C. THE BUGS I SHIPPED, and the rules that would have caught them

Four, in two families. All four were mine and three of them silently produced confident wrong
numbers rather than errors.

## C1. Cross-artifact lookups without a known-good test case (twice)
* **Pair-order:** I assumed `pins.json` stores `(x, y)`. It stores `(y, x)`. My first lookup
  returned NOT FOUND for *every* variable — including x24601/x2081, whose answers four agents
  already knew.
* **Reduction frame:** wire values are **unreduced ~89-digit integers**; I searched for the
  reduced/shifted coordinate and reported "**0 of 256** ladder points on wires", which was a false
  negative that would have over-killed my own model for the wrong reason.

> **Rule: never join two artifacts without first running the join on a pair whose answer is already
> known.** `x24601 -> 72` and `x2081 -> 235` are this lab's canonical test; a join that fails them
> is broken regardless of how plausible its output looks.

## C2. Filters that silently drop data
* **Regex:** `relax.py` used `re.fullmatch(r'\(?x\d+\s*\*?.*')`, which fails on expressions
  starting `((` — silently dropping atom **7887** (`x24267`, cost 4). That atom turned out to be
  half of the best floor I later found. My "exhaustive" pair scan had also only ranked the top 25.
* **Index namespaces:** I cross-quoted atom IDs from `NOTEBOOK.md` §Session 10 (22229, 35758…)
  against `E`'s indexing (3130, 7251…). Different numbering; the comparison was meaningless. This
  has now caught four agents.

> **Rule: when a filter selects a subset, print how many it dropped and reconcile against the
> exhaustive count. When quoting an index, name the artifact that defines it.**

## C3. The meta-lesson, which is the one that actually cost me
My `fold(k) != T` "validation" passed and meant almost nothing — **almost any wrong model also
yields `fold != T`.** A validation that a wrong model would also pass is not a validation.
The test that finally worked was checking predicted **wire contents** against the one verified
object. **Validate a model by what it predicts is PRESENT, not by what it predicts is absent.**

---
# §D. Files
Incidence work (§A, durable): `footprints.py rank.py reach2.py cancel.py obstruct.py price.py
defect.py cfgscan.py rescore.py crosscheck.py` -> `runs/{footprints1,rank,reach2,cancel,obstruct,
price,defect,cfgscan,rescore,crosscheck}.json`
Validation (§A7, §B1): `validate_A.py validate_A2.py validate_A3.py depth.py` -> `runs/validate_A*.json`
Withdrawn (§B): `solve2.py collapse.py tradeoff.py relax.py realize.py supports.py`
Model/derivation (§B3): `model.py group.py ladder.py order.py fastgrp.py ladder.json points_short.json`
Solver benchmarks (§A8): `z3enc.py encode.py witness.py sibling.py bench.py cnfbench.py
opt3_cpsat.py opt2_uf.py search_lw.py bsgs.py` -> `encodings/`, `runs/{bench,cnf,cnfsize,cpsat}.json`
Narrative: `LOG.md`. Previous revision of this file: `runs/RESUME_R_prev.md`.

# §E. Standing rules honoured
- Every "nothing can move X" statement carries its knob set and selector configuration
  (`LOG.md` §7, §9, and the SCOPE block in §9).
- Nothing I produced exceeds 4,300 digits, so `checker.py` was valid for everything I verified.
- No generator forensics and no named-curve framing: `B`, the ladder and `N` are measured from the
  decoded law by the four scripts in §D and are stated in §B3 as identities about a point set.
- Reads outside my directory were confined to read-only `agentF_work` imports plus the shared
  `checker.py` / `best/`. No git commands were run at any point.
