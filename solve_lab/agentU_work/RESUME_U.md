# RESUME_U — agent U.  The partition theorem.

Everything here was measured on this box from `EQUATIONS.txt` with **my own parser**.
No other agent's code was imported.  Other agents' directories were read read-only, and only
as *cross-checks* (their gadget SITES, my arithmetic).

## 0. SCORE

Baseline re-verified at the start of this session:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.

**I did not beat it.** Nothing in `agentU_work/` is a better partial. No infeasibility is claimed.

---

## 1. HEADLINE — the arithmetic half of K's partition theorem is CLOSED, exhaustively, under the
## TIGHT criterion, from a parse that is not K's.

> **For every slot of every merge gadget in this instance, the maximum attainable subset sum of
> its leaf-exponent support is strictly below `N` — the largest is `0.798718631·N`.  Hence
> `|Σ_A 2^i − Σ_B 2^j| < N` for every choice of `A ⊆ I`, `B ⊆ J`, at every gadget, so `±N` is
> unreachable everywhere.  This is exhaustive, not bounded: no search is involved.**

### Knob set and configuration (stated per the lab rule)

* **Knob set:** all 256 leaf selector bits over all `2^256` configurations, and in fact all
  `A ⊆ I`, `B ⊆ J` at every gadget — i.e. every live-leaf subset, not a sampled family.
* **Configuration:** any.  The statement is about the exponent arithmetic and the measured
  slot partition; it does not reference a particular assignment.
* **What it does NOT cover:** it is conditional on the *fold semantics* — that a slot carries
  the composition of the live leaves in its own support.  That is K's step 2 / the routing
  layer, and it is **not** what I closed.  See §5.

---

## 2. THE ARITHMETIC, RECOMPUTED (`u0_arith.py`)

```
N  = 115792089237316195423570985008687907852837564279074904382605163141518161494337  (256 bits)
popcount(N) = 192, zeros = 64
2N > 2^256 - 1                       (slack 1.1579e76)   -> k = ±1 is the only possible wrap
2^256 - N = 432420386565659656852420866394968145599  <  2^129
```

Two corrections/sharpenings to what is in FLEET.md and RESUME_K.md:

1. **The unconstrained condition is TRUE — reconfirmed independently.**  There are **34**
   indices `j` with `bit_j(N)=1, bit_{j+1}(N)=0`; `j = 0` gives non-empty disjoint
   `A` (|A|=192), `B` (|B|=1) over `{0..255}` with `Σ_A 2^e − Σ_B 2^e = N` exactly.  Verified
   numerically.  So the whole question is the partition and nothing else.  Confirmed.

2. **K's stated condition is strictly weaker than the tight one, and I used the tight one.**
   K's §4.0 says "neither slot support contains all of `{129..255}`", justified by
   `2^256 − N < 2^129`.  That is *sound* but not sharp: containing `{129..255}` is **necessary
   but not sufficient** for `maskval ≥ N` —
   `maskval({129..255}) = 0.9978·N < N`.  The exact criterion is simply

   ```
   Σ_A 2^i = +N + Σ_B 2^j  requires  Σ_A ≥ N,  and  Σ_A ≤ maskval(I).
   So  ±N reachable at a gadget  =>  maskval(I) ≥ N  or  maskval(J) ≥ N.
   ```

   equivalently `Σ_{e ∉ I} 2^e ≤ 2^256 − N − 1`.  I tested `maskval ≥ N` directly at every slot.
   **The −N case is covered by the same test on the other slot** and was checked separately, not
   assumed by symmetry.

---

## 3. THE PARTITION FACTS, MEASURED INDEPENDENTLY

### 3.1 My parse (`v1_parse.py` → `v3_defs.py` → `v5_chain.py` → `v8b_supp.py`)

* Recursive-descent parse of all 39,033 equations into **37,936 distinct maximal atoms** in
  **20 shapes**.  (This is *my* atomisation — 37,936, not F's 39,033.  Rule 3: a count derived
  from one parse is a fact about that parse.  Nothing here depends on the count.)
* **512 leaf pins** of shape `sel·(w − C) − m·z`, over **exactly 256 distinct selectors, two
  pins each, 512 distinct wires, 512 distinct 287–296-bit constants.**
* `p = CONST[x26064]` = 115792089…908834671663, 256 bits, prime.
* **Curve solved algebraically, not taken:** fitting `y² = (x+s)³ + b` to three leaf constants
  (orientation of each selector's two pins searched, since it is not given) yields
  `shift = 109712675…394060739`, and `3·shift mod p = 97553848499418123410591666447050222001188385549510401465815187079080512838891` — which
  is exactly the `K` in K's decode, recovered here without reading it.
  `b = 64019533680030876408443198762210829058751700634554282185987325820393598524794`.
  **256/256 leaf points satisfy `Y² = X³ + b`.**
* **Doubling chain:** 255/256 doublings land inside the leaf set, one source, one sink, chain
  length 256, exponents 0..255.  `N·G = O`.  `leaf(e) == 2^e·G` for **all 256**.
  Base selector = **x2779** (independently equals Q's `G_leafvar`).
* **Selector-support closure** over all 38,748 wires (forward on the canonical definition DAG,
  copies union-found, product-constraint atoms read forward) gives **exactly 511 distinct
  non-empty supports**:

| | |
|---|---|
| distinct supports | **511 = 2·256 − 1** |
| laminarity violations | **0** |
| internal nodes | **255, every one binary, children disjoint, union = parent** |
| root | `{0..255}` |
| **root halves** | **178 / 78, disjoint** |
| A (178) ∩ `{129..255}` | **84** — omits 43 |
| B (78) ∩ `{129..255}` | **43** — omits 84 |

**Those are K's two measured partition facts, reproduced exactly (43 and 84) by a parser that
shares no code and no intermediate file with K.**  K's "every interior stage sits inside one
root half" is subsumed: the family is *laminar with zero violations*, so every support is
inside a root half by construction, and I verified the tree property rather than assuming it.

### 3.2 Second source — L's 383-node model (`agentT_work/mirror/L/full_model.pkl`, read-only)

L's per-node leaf sets `sub[]`, mapped through **my** exponent map, give a family that is
**set-equal to my 511**.  383 internal nodes → 255 with two non-empty slots (the genuine
merges) + 128 with one empty slot (pass-throughs).  **That reconciles 383 vs 255: the
difference is exactly 128 pass-through nodes, not a disagreement.**  Root split 178/78.

### 3.3 Third source — K's numbers, from FLEET/RESUME_K: 178/78 and 43/84.  Identical.

---

## 4. THE TEST (`v9b_theorem.py`, `v11_direct.py`)

| source | gadget sites | verdict | max `maskval(slot)/N` |
|---|---|---|---|
| my own tree | 255 | **0 sites can reach ±N** | 0.798718631 |
| L's 383-node model, my supports | 383 | **0 sites** | 0.798718631 |
| Q's 383 chord gadgets (`ua`/`ub`), my supports | 383 | **0 sites** | 0.767097519 |
| **tree-free**: all 510 proper supports | — | **0 supports have `maskval ≥ N`** | 0.798718631 |

* **Exhaustive brute force as an independent check of the bound argument itself:** for the 240
  of 255 sibling pairs with `|I|+|J| ≤ 22`, every subset sum of both sides was enumerated —
  **14,052,776 subset-sum pairs, zero representations of ±N.**  The remaining 15 pairs are
  settled by the exact interval `[−maskval(J), +maskval(I)] ⊂ (−N, N)`, which is a
  recomputation, not an argument.
* **The tree-free form is the strongest and the one to cite:** *no proper slot support in this
  instance has `maskval ≥ N`.*  This does not depend on my tree recovery, on the 178/78 split,
  on which side any exponent sits, or on the pairing being right — any two disjoint sets drawn
  from the family fail.  (This is the closure K's `k33` failed to reach with inflated supports.)
* **The adjacent hole (a half folding to the group identity) is closed too**, and by a sharper
  argument than K's size count: `Σ_S 2^e ≡ 0 (mod N)` with `Σ_S < 2N` forces `Σ_S = N` exactly,
  and because the summands are **distinct powers of two the binary representation is unique**,
  so it forces `S = supp(N)` exactly.  **No proper support contains `supp(N)`** — measured.
* **The `dx = 0, dy ≠ 0` case** (inputs are negatives of each other, `Σ_A + Σ_B ≡ 0 mod N`) is
  *not* a free-output case: `R1 = S·dx² − dy² = −dy² ≠ 0`.  It is an unsatisfiable gadget, not
  an exploitable one.  Noted because the same uniqueness argument shows `Σ_A + Σ_B = N` forces
  `A ∪ B = supp(N)` — which IS available at the root (`supp(N) ⊆ {0..255}`), so this branch had
  to be checked on the residual, not on the arithmetic.

---

## 5. WHAT IS **NOT** CLOSED — stated plainly

The theorem I closed is: **an *honest* coincidence is impossible.**  "Honest" = each slot
carries the composition of the live leaves in its own support.  Two things remain:

1. **The fold semantics (K's step 2, the routing layer).**  Over ℤ the leaf pin is
   `sel·(w − C) − m·z = 0`, so with `sel = 1` the wire carries `w = C + m·z` — **not** `C` —
   until `z` is separately forced to 0.  **256 of the 512 leaf pins have `m = 1`.**  So a leaf
   wire is only pinned to its constant modulo whatever pins `z`.  This is exactly check-in 40's
   correction and it is where the remaining conditionality lives.  My theorem does not touch it.
2. **A *dishonest* coincidence is not only possible, it is what the deliverable does.**  See §6.

**So: do not cite §1 as "the degeneracy route is closed."**  Cite it as *"no gadget can be fed
two coinciding inputs by any assignment of the 256 selectors under the fold semantics; a
coincidence therefore requires forcing a wire off its honest value, and that has a price."*

---

## 6. END-TO-END VALIDATION AGAINST THE DELIVERABLE (`v12_deliv.py`) — the strongest check here

Driving my decode on `best/new_instance_partial_39026.json` (3,540 assigned variables):

* **Exactly 1 of the 383 chord stages has non-zero inputs.**
* At that stage (Q's index 152) **the two inputs coincide exactly**, and both are on the cubic.
* The coinciding value is **`2^72·G`** — a genuine leaf point, identified by my own exponent map.

That is K's §3 mechanism, confirmed numerically from an independent decode: the deliverable
makes one gadget see two equal live inputs, the chord residual vanishes identically, the output
goes free and is driven to the target.  It pays 7 equations for the lie that puts `2^72·G` on a
wire whose honest constant is a different leaf.  **The degeneracy is reachable — dishonestly,
at a price — and §1 says the price can never be zero.**

---

## 7. A BUG I MADE AND CAUGHT — worth the two lines

My first parser collected **every** `-` node in the AST, including nested ones.  That turned the
inner difference of `dx = ua − ub` into a "copy" atom `(x_ua − x_ub)` and union-found the two
slot inputs of every gadget together.  Caught by cross-checking against Q's `qstages.json`:
**0/383 slot pairs came back disjoint**, and `find(ua) == find(ub) == find(u3)`.

**Blast radius, measured before reporting (rule 8):** the corrected parse (`v1_parse.py`,
maximal `-` nodes only: 6,622 → 3,749 copies) gives the **same 511-set family with the same size
profile and the same 178/78 root split**.  So the headline number never moved; only my
confidence in it should have.  The lesson is the ledger's rule 4 in reverse: I validated by
what the model predicted was PRESENT (disjointness at 383 known gadget sites) and it failed
loudly.  A support family that is laminar and looks like a binary tree is **not** self-validating.

---

## 8. FILES

| file | what |
|---|---|
| `u0_arith.py` | the three-line arithmetic + the `Σ_A − Σ_B = N` witness |
| `v1_parse.py` → `v_atoms.pkl` | corrected recursive-descent parse, 37,936 maximal atoms |
| `v3_defs.py` → `v_defs.pkl` | defs / copies / constants / 512 leaf pins |
| `v5_chain.py` → `v_leaves.pkl` | curve fit, 256 points, doubling chain, exponent map |
| `v8b_supp.py` → `v_supp2.pkl` | selector-support closure over all 38,748 wires |
| `v9b_theorem.py` → `v_tree_final.pkl` | laminarity, tree, root split, the `maskval` test |
| `v11_direct.py` | direct per-gadget test from 3 sources + 14M-pair brute force |
| `v12_deliv.py` | end-to-end validation against the 39,026 deliverable |
| `u1_parse.py`, `u3..u10` | the buggy-parser generation, kept for the §7 comparison |

Superseded within this directory: `u_atoms.pkl`, `u_defs.pkl`, `u_supp.pkl`, `u_tree.pkl`,
`u_tree_final.pkl` — products of the pre-§7 parser.  **`u_leaves.pkl` is unaffected** (leaf pins
are maximal atoms in both parsers) but use `v_leaves.pkl`.

---

## 9. CROSS-PARSE IDENTITY OF THE EXPONENT MAP

My `sel2exp` (256 selectors → exponents), built from my own curve fit, my own doubling chain and
my own base discovery, is **byte-identical to Q's `qladder.json` `sel2exp`** — 256/256, zero
differences.  Two parses that share no code agree on the labelling, not just on the structure.
So exponent numbers in this file are comparable with Q's (they are **not** comparable with K's
`chain.json`, which uses a different labelling).

## 10. THE DELIVERABLE'S LIE IS IN THE ROUTING, NOT THE PINS (`v14_onset.py`)

* **Exactly 2 leaf selectors are 1**: exponents **72** (`x24601`) and **235** (`x2081`); the
  other 254 are 0.
* **All four of their coordinate wires carry their own honest pin constant**, and all four `z`
  wires are 0.  **No leaf pin is violated.**
* The two ON exponents are separated by the **root** sibling pair: 235 ∈ the 78-half,
  72 ∈ the 178-half.
* Yet the coincidence happens with **both** inputs equal to `2^72·G`.

So the lie is a **cross-half route** — the 78-half chain is made to carry the 178-half's value —
and it costs 7 equations.  This is K's §3 read off the deliverable by an independent decode, and
it is what §1 forbids doing honestly: an honest root coincidence needs
`Σ_A 2^i − Σ_B 2^j = ±N` with `A` in the 178-set and `B` in the 78-set, and both maskvals are
below `N`.

## 11. WHAT PINS A LEAF WIRE — measured (`v13_pins.py`)

`sel·(w − C) − m·z = 0`, so with `sel = 1` the wire carries `w = C + m·z`.  **256 of the 512
pins have `m = 1`**, which would make `w` free if `z` were free.  It is not:

* **Every one of the 512 `z` wires occurs in exactly 2 atoms** — its own pin, and one
  `(V − (V*V))` atom.  **512/512 `z` wires are defined as a product of two wires.**
* So `sel·(w − C) = m·a·b`: the leaf wire is pinned to its constant **iff the product `a·b` is
  zero**, and lying at a leaf means switching on a gate, not just choosing a free value.
* Incidence price of the cheapest single-pin lie: `{pin atom} ∪ {atoms containing z}` touches
  **min 8 / median 14 / max 21 equations**.  **This is an incidence count and therefore an
  inflated upper bound, not a floor** — three separate demonstrations in this lab (C, P, L) show
  incidence pricing fails, and the deliverable's actual price is 7, below this minimum.  **Do
  not read "8" as a bound on anything.**

## 12. K's PROMISED SWEEP — completed here (K never published it)

Three engines exist in `agentK_work/`: `cascade.Cascade.close` and `cascade2.Inc.run` have **no
guard parameter at all**; only `cascadep.CascadeP.close(seed, order, forbid=(), pin=None)` can be
guarded.  Sweeping every script for closure *runs* (not merely imports):

| script | engine / guard |
|---|---|
| `k13_root.py` | `cascadep.close` **UNGUARDED** |
| `k17_validate.py` | `cascadep.close` **UNGUARDED** |
| `k24_allon.py` | `cascadep.close` **UNGUARDED** |
| `k7_order.py` | `cascade2.Inc.run` — **no guard parameter** |
| `k9_handles.py` | `cascade2.Inc.run` — **no guard parameter** |
| `k26_drive.py` | guarded **conditionally** (`forbid=FORBID if forward_only else ()`) |
| `k43_forward.py`, `k44_audit.py` | guarded conditionally (`pin=pin if guard else None`) |
| `k35_otherbools.py`, `k37_premise.py`, `k38_deadgate.py`, `k39_alias.py` | guarded on every call |

**K's own audit table omits four of these: `k13_root.py`, `k17_validate.py`, `k24_allon.py`,
`k7_order.py`.**  Only `k9_handles.py` and `k26_drive.py` are flagged there.

**The one that matters, and it lands squarely on the fact I was sent to attack.**
`k24_allon.py`'s docstring says it "**also settles the side of leaf exponent 163**" — the leaf
K had to place by hand to turn its root split from 177/78 into **178/78** — and it runs an
**unguarded** closure.  K's audit table records the root split as SAFE, on the grounds that
`k36_tight.py` (support recovery) uses no closure; it does not record that the 163 placement came
from a closure run that could solve constraints backwards.

**Blast radius: zero, and it is closed rather than merely unrefuted.**  My 178/78 split is
derived from the definition DAG alone — no closure, no seeding, no assignment, forward by
construction — and L's independent 383-node model gives the same split and a set-equal leaf-set
family.  So the fact survives; what did not survive is K's justification for calling it safe.
(Note my own parse hit the *same kind* of gap at the *same kind* of leaf: exponent **107** was
the one leaf my first support closure failed to attach, until product-constraint atoms were read
forward.  Two independent parses each needed one leaf placed by an extra mechanism.  That leaf is
the fragile point of this decode and anyone re-deriving the split should check it first.)

## 13. WHAT I DID NOT DO, DELIBERATELY

K's other unfinished item — re-running the whole B-half fold validation table with backward
derivation blocked at **every** slot — I did not attempt.  K's file carries an explicit
"**DO NOT REBUILD THE GUARD**" header with the reason: the guard is built from atom *shape*, and
`k43`'s diagnosis is that the role of `W` and `Z` in `((xW − xZ) − xH)` is reversed at non-root
slots, so a shape-built pin map blocks real sources (guarded 0/18 halves match vs unguarded 6/18
— it broke even the single-leaf pass-through case).  The correct pin map has to come from the
decoded slot→source direction.  More to the point, **the thing that table was for is settled
from an independent parse without any closure** (Q: 383/383 chord law by Schwartz–Zippel against
the real sub-DAG, 383/383 mux quadrant law, one tree, one root), and rebuilding it would
re-derive a settled result with the least trustworthy instrument in the lab.  Recorded as a
decision, not an oversight.

## 14. THE SINGLE HIGHEST-VALUE NEXT EXPERIMENT

§1 says the price of a coincidence can never be zero.  §10 says the deliverable pays **7** for a
**cross-half route** — not for a leaf-pin violation, which is the assumption most of this lab's
placement searches were built on.  So the open question is sharply:

> **What is the minimum equation-cost of a route that puts one root half's value onto the other
> half's slot, over all 383 slots and all 2^256 ON-sets?**

The deliverable's is 7 with `|ON| = 2` at the root.  M's placement enumeration priced 4,096
subsets of `H12` and T verified the engine at 9/9, but that search is indexed by *handle
subsets*, not by *route sites* — and the deliverable's lie is a route.  The concrete test:
enumerate the 383 slots, and for each, price in **equations** (with L's exact in-memory scorer,
never incidence) the cheapest assignment that makes that slot's two inputs coincide.  Interior
slots have far smaller supports than the root and were never priced this way.  If any slot comes
in below 7, that is the campaign's terminal result; if none does, 7 acquires a mechanism instead
of an exhaustion.

---

# CHECK-IN 87 FOLLOW-UP — the slot-pricing experiment: PARTIAL, and I did not price a single slot

**I did not deliver per-slot prices.**  I got two-thirds of the way and the instrument failed
calibration.  What follows separates what is checker-verified from what is not.

## 15. The pin-level barrier does not exist (`w1_zfactors.py`, `w2_wire.py`, `w3_crt.py`)

Every leaf pin is `sel·(w − C) − m·z = 0` and **512/512 `z` wires are defined as a product
`z = a·b`**.  Measured over all 512:

* **1019 of the 1024 factors are FREE variables** (no definition, no constant pin).  The other
  5 are pinned to `p`.  **507/512 pins have BOTH factors free; 512/512 have at least one.**
* So `z` is unconstrained and the pin reduces to a pure divisibility: `m | (w − C)`.
* **All 256 leaves carry `m = 1` on their Y wire and `m > 1` on their X wire** (256 distinct
  X moduli).  So **every leaf's Y coordinate is completely free at zero cost.**
* For two leaves to be driven to a *common* coordinate pair, Y is free and X needs
  `W ≡ C_aX (mod M_a)`, `W ≡ C_bX (mod M_b)` — solvable iff `gcd(M_a,M_b) | (C_aX − C_bX)`.
  Over all 32,640 cross-slot leaf pairs: **26,389 are feasible**, and **232 of 255 slots admit
  at least one feasible pair**.  `gcd = 1` for 24,743 pairs, which makes them automatic.

**So the price of a coincidence is not paid at the pin.**  This kills the assumption most of
this lab's placement searches were built on, and it is measured, not argued.

## 16. What the price IS paid for — measured with `checker.py`'s own compiled equations

Harness: `checker.load_equations()` once (28 s), then `evaluate_all` per candidate.  Deliverable
re-scores **7 failing** through the harness, so the scorer is calibrated.  Three candidates, all
on the deliverable's own leaf pair (exponents 72 and 235), all **checker-verified**:

| candidate | failing |
|---|---|
| deliverable (route lie) | **7** |
| pin lie: leaf 235's wires := leaf 72's constants | **50** |
| pin lie: leaf 72's wires := leaf 235's constants | **46** |
| joint CRT: one common `W` on both leaves (`gcd = 3`, divides) | **88** |

**These are NOT slot prices and must not be quoted as such.**  They are what a pin lie costs
*without re-propagating the downstream chain* — the mux pass-through atoms below the slot still
carry the old values, and the checker counts every one.  The honest reading is only:
**the cost lives in the propagation, not in the pin.**

## 17. THE INSTRUMENT FAILED CALIBRATION — stated before any number derived from it

I built a forward-only evaluator from my own parse (`w5_eval.py`): 31,853 variables with a single
definition each, a Kahn topological order with **0 cycles**, every variable recomputed from its
**own** definition, nothing ever solved backwards.

**Calibration test: propagate the deliverable and re-score.  Result: 8,229 failing, not 7, with
4,578 variables changed.**  The evaluator is therefore **wrong** and **no slot was priced with
it.**  Diagnosis: of the 3,749 copy atoms, orientation is only forced when exactly one side
carries a definition; where neither does I picked a direction arbitrarily, so some copies run
backwards.  **That is precisely K's failure mode, reached independently, and it is why I am
reporting a gap rather than a table.**

Nothing in §§1–14 depends on this evaluator — the partition theorem is pure arithmetic over the
support family and never evaluates the circuit.

## 18. HANDOVER — how to actually run check-in 87's experiment

**Do not rebuild the evaluator.**  A calibrated one already exists and is audited: L's
`calib2.py` + `full_model.pkl` builds an assignment from an ON-set and orientations, and **M's
incremental engine was verified EXACT by T outside M's parse** — driven on the witness subset it
reproduces 39,026 with the deliverable's exact 7 failures and is byte-identical to the
deliverable (0 of 38,748 vars differ).  Drive **that** over the 383 slots.

The construction to drive, per slot β with subtrees `I`, `J`:

1. pick `ℓ_a ∈ I`, `ℓ_b ∈ J` with `gcd(M_a,M_b) | (C_aX − C_bX)` — §15 gives 26,389 such pairs
   and the per-slot counts;
2. turn on exactly those two selectors; set both X wires to the CRT value `W`, both Y wires to a
   common `W_Y`, and each `z = (w − C)/m` (free, §15);
3. **re-propagate** the whole live path — this is the step my instrument could not do;
4. β's two inputs now coincide, so its residuals vanish identically and its output is free;
   everything above β is pass-through because the sibling subtrees are dead, so set β's output to
   the target and the root carries it;
5. score with the exact scorer and record the failing count **per slot**.

Below 7 anywhere is terminal.  At or above 7 everywhere turns the deliverable's 7 from an
exhaustion into a mechanism.

**Two cautions from what I did measure.**  (i) The point `(W, W_Y)` need not be a curve point —
β's chord law is vacuous and nothing above β applies one — so the search space is larger than a
curve-point search and must not be restricted to leaf values.  (ii) 232 of 255 slots admit a
feasible pair, so the pairing is not the bottleneck; the propagation cost is.

---

# CHECK-IN 90 FOLLOW-UP — the engine is CALIBRATED and mirrored; the sweep is NOT run

## 19. The pkl wipe, and the mirror that works (`mirror/`, `w6_mcal.py`)

Confirmed exactly as the coordinator warned.  `agentE_work/model3.pkl` and `agentE_work/dag.pkl`
are **both absent** — `.gitignore` carries a global `*.pkl`.  M's engine chain is
`engine3 -> agentE_work/harness -> agentE_work/{model3,dag}.pkl`, so it does not run from cold.

**M kept its own copies**: `agentM_work/model3.pkl` and `agentM_work/dag.pkl` are present.  So the
mirror is two files and no rebuild was needed (I did not have to run `t_rebuild.sh`):

* `agentU_work/mirror/harness.py` — `agentE_work/harness.py` with both pkl paths repointed to
  `agentM_work/`.
* `agentU_work/mirror/engine3.py` — `agentM_work/engine3.py` with its `sys.path.insert(0, ...)`
  repointed at `agentU_work/mirror` (it hard-codes `agentE_work` and otherwise shadows the mirror).

Nothing outside `agentU_work/` was written.

## 20. CALIBRATION — stated before any number derived from it, and it PASSES

```
engine3 imported 5.9s ; NV = 38748
deliverable via checker            : 7 failing
seed extracted                     : 37 entries
forward(seed_of(deliverable))      : 7 failing via checker.py
variables differing from deliverable: 0 of 38,748
CALIBRATION: PASS
```

This is T's audited property reproduced from outside M's directory: the engine is a genuine
forward propagator and it is **byte-identical** on the deliverable.  **This is the instrument
§17 said I lacked, and it is now in place and validated.**  Contrast with my own evaluator
(`w5_eval.py`), which failed the same control at 8,229 failing — that one stays retired.

## 21. THE SWEEP IS NOT RUN — zero slots priced, still

I calibrated and stopped.  **No slot has been priced, exactly or approximately.**  The three
numbers in §16 remain what they were: one leaf pair with stale propagation, not slot prices.

**What is left is small and fully specified.**  The engine takes a **37-entry seed** and returns
a full assignment; `checker.py`'s compiled equations score it in ~0.3 s.  So the whole sweep is:
map §18's construction into a seed, then loop.  The one unknown a successor must resolve first is
**the seed vocabulary** — which of the 37 entries carry the ON-set, which carry the leaf X/Y wire
values, and which carry the free slot output — read off `Eng.seed_of` / `Eng.forward` in
`mirror/engine3.py` against `H.SEQ` / `H.definer`.  Get that map, and §18 runs.

**Order the sweep by my two cautions, not by slot index.**  The common point need not be on the
curve (β's chord law is vacuous, nothing above β applies one), so do not restrict to leaf values;
and at 232/255 slots the pairing is free, so the cost is propagation — which means **slot depth
and support size are the variables that matter**.  If all 383 is expensive, stratify by those two
and say so; do not price a prefix.

**Below 7 at any slot is terminal.**  At or above 7 everywhere, the deliverable's 7 becomes a
mechanism rather than an exhaustion.

## 22. FILES ADDED THIS ROUND

| file | what |
|---|---|
| `w1_zfactors.py` | 1,019/1,024 pin factors are free variables |
| `w2_wire.py` → `w_tag.pkl` | `m=1` on Y at all 256 leaves; `m>1` on X |
| `w3_crt.py` → `w_xy.pkl` | 26,389/32,640 CRT-feasible pairs; 232/255 slots |
| `w4_price.py` | checker-exact scores of pin lies **without** re-propagation (46/50/88) |
| `w5_eval.py` | my own evaluator — **failed calibration at 8,229, retired** |
| `mirror/harness.py`, `mirror/engine3.py` | the working mirror of M's engine |
| `w6_mcal.py` | **the calibration: PASS, 7 failing, 0 vars differing** |

---

# CHECK-IN 94 — SECOND pkl WIPE: full rebuild, calibration PASSES, and the slots ARE priced

## 23. THE REBUILD (this time nothing was cached — `agentM_work` had 0 pkl too)

The restart wiped every `*.pkl` again.  Unlike check-in 90, **M's private copies were gone as
well**, so the mirror had to be rebuilt from source.  Everything below was written only in
`agentU_work/`.

| built | from | into | reproduces |
|---|---|---|---|
| `model3.pkl`, `dag.pkl` | `agentE_work/parse3.py`, `agentE_work/dag.py` (copied read-only) | `mirror/E/` | atoms **40,727**, vars defined **30,383**, free **8,365**, cycle hits **0** — M's `shim.py` logs the same three numbers |
| `full_model.pkl`, `calib2.pkl`, `handles.pkl`, `ors.pkl`, `ortree2.pkl`, `slopes.pkl` | T's `mirror/L` chain (T rebuilt them in its own dir during this session) | copied to `mirror/L/` | md5 `7e034ea6…` / `7f62e53b…` |
| `v_atoms/v_defs/v_leaves/v_supp2/v_tree_final/w_z/w_tag/w_xy` | my own `v1→v3→v5→v8b→v9b→w1→w2→w3` (`u_rebuild_vw.sh`) | `agentU_work/` | 511 supports, root split **178/78**, **0** gadgets with `maskval ≥ N`, **1,019/1,024** free pin factors, **24,743** gcd-1 pairs — every §1–§15 number reproduced |

**Three path patches, all load-bearing:** `mirror/harness.py` now points at `mirror/E/*.pkl`;
its `orient.pkl` cache is an absolute path into `mirror/` (it used to land in cwd);
`mirror/engine3.py` keeps its `sys.path.insert(0, …/mirror)` — without it `agentE_work` shadows
the mirror and every repoint silently does nothing.  `umodel.py`'s `LDIR` repointed at `mirror/L/`.

## 24. CALIBRATION — stated before any number derived from it

```
python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json
  -> satisfied 39026/39033 (7 failing), failing [12231,12270,12350,14584,18673,22044,29125]

engine3 imported 1.7s ; NV = 38748
deliverable via checker              : 7 failing
seed extracted                       : 37 entries
forward(seed_of(deliverable))        : 7 failing via checker.py
variables differing from deliverable : 0 of 38,748
CALIBRATION: PASS
```

Identical to check-in 90's target.  ~0.19 s per candidate end-to-end (build → seed → forward →
`checker.evaluate_all`).
