# LOG_R — agent R: automated reasoning against the REDUCED problem

Angle: SAT / SMT / CP / MIP against the decoded instance, not against the 39,033 raw equations.
Agent C's earlier abandonment was correct *for the undecoded instance*; this log re-opens the
question against agent F's decode and answers it with measurements.

## 1. Baseline re-verified (mine, from cold)
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`, 26 s.

## 2. Before encoding anything: what exactly is the reduced problem?
The brief said to encode "96 stages, 256 boolean selectors, one uniform degree-3 law, one target
at the root". Encoding that literally would have been wrong, because two of those four are not
binding. Measured here (`model.py`, `group.py`, `ladder.py`, `order.py`):

1. **The offset disappears.** Substituting `X = x + K/3` (3 is invertible mod p) turns
   `out_x = l^2 - a_x - b_x - K` into `out_X = l^2 - a_X - b_X`. Checked against `chordK` on 200
   random leaf pairs: identical.
2. **All 248 forced pin points and the target satisfy one relation `Y^2 = X^3 + B`**, with
   `B = 64019533680030876408443198762210829058751700634554282185987325820393598524794`,
   fitted from 2 points and then verified on the other 246 plus the target: **0 exceptions**.
3. **The fold law is therefore an abelian group law** — commutativity and associativity checked
   on 200 random triples each. **So the fold of a leaf subset does not depend on the tree shape.**
   The 96 stages are *not* 96 independent constraints; they are one commutative sum. F's
   outstanding decode work (56 undecoded slot pairs, 24 leaf-adjacent stages) is **not needed**
   to state or to attack the problem.
4. **The 256 leaves are one doubling chain.** `ladder.py`: 239 of the 248 forced points have their
   double also in the set; the set splits into 9 doubling-closed pieces (111, 61, 28, 12, 11, 9,
   8, 5, 3 = 248), and each piece's last point, doubled twice, lands exactly on the next piece's
   first point — 8 splices, one unknown point per gap, one head. So the 256 leaves are
   `L_i = 2^i * L_0`, i = 0..255, and the 8 points `sweep_ii` failed to force are recovered as
   the gap points. Checked: `L_i == 2^i * L_0` for sampled i.
5. **The group is cyclic of 256-bit prime order.** `order.py`, Cornacchia on `4p = L^2 + 27M^2`,
   then the unique candidate that annihilates the base point:
   `N = 115792089237316195423570985008687907852837564279074904382605163141518161494337`,
   prime (`sympy.isprime`), 256 bits. Prime order ⇒ no Pohlig–Hellman decomposition exists.

**Conclusion (the encoding target):** the reduced problem is
> find the 256-bit scalar `k` with `k * L_0 = T`, bit *i* of `k` being the selector of leaf `L_i`.

### Validation demanded by RESUME_F — PASSED
Deliverable ON-set `{24601, 2081}` -> ladder indices **{72, 235}** -> `k = 2^72 + 2^235`, and
`fold(k) != T`. That is exactly the prediction F said a correct evaluator must produce (it must
*not* reproduce the target). The deliverable pins the root wires to `T` while its own leaves fold
elsewhere; that gap is the 7 failing equations.

## 3. Encodings built (all on disk, re-runnable without re-encoding)
`z3enc.py` (z3py builder, Int and BV), `encode.py` (raw SMT-LIB writer), `witness.py`
(independent Python mirror that computes the intended solution and asserts every encoded
relation — used to prove each encoding is satisfiable before any solver is blamed),
`opt3_cpsat.py` (OR-Tools CP-SAT), `opt2_uf.py` (uninterpreted-law variant),
`sibling.py` (scaled-down instances of the *identical shape* over an m-bit prime).

Per stage the encoding is division-free (no modular inverse variable):
```
d  = bx - ax          n  = by - ay            (guard d != 0)
(sx + ax + bx) * d^2 == n^2                   (chord, x)
(sy + ay) * d        == n * (ax - sx)         (chord, y)
acc' = acc + b_i * (S - acc)                  (mux)
```
7 modular multiplications + 1 boolean per stage; accumulator seeded at `L_0` so the identity
element never appears. `witness.py` confirms the intended solution satisfies all of it at
m = 8, 12, 16, 20.

## 4. MEASURED RESULTS  (4 cores, load 13–17 from the rest of the fleet)
See `runs/bench.json` (z3), `runs/cpsat.json` (CP-SAT), `runs/opt2.log`.

### 4a. SMT — z3 5.0.0, 300 s limit  (`runs/bench.json`)

| instance | logic | selectors | result | time |
|---|---|---|---|---|
| m=8 (p=181, 7 stages) | QF_BV | **all pinned to the known solution** | sat | **119.1 s** |
| m=8 | QF_BV | free (128 possibilities) | **timeout** | 300 s |
| m=8 | QF_NIA | all pinned | **timeout** | 300 s |
| m=8 | QF_NIA | free | **timeout** | 300 s |

The pinned rows are the important ones: with every selector fixed the instance is a
straight-line computation with a unique solution, and z3 still needs 119 s (QF_BV) or fails
outright (QF_NIA). Both encodings are *proved satisfiable* first by `witness.py`, so these are
solver failures, not encoding failures.

### 4b. CNF — z3 bit-blast -> CaDiCaL 1.9.5 / Glucose 4.2 (`runs/cnf.json`, `encodings/*.cnf`)

m=8 pinned: 57,299 vars / 388,439 clauses.
* **CaDiCaL 1.9.5**: sat, **107.6 s**, 165,725 conflicts, 354,424 decisions, 34,556,815 propagations.
* **Glucose 4.2**: sat, **276.7 s**, 353,100 conflicts, 535,750 decisions, 122,898,960 propagations.
* Kissat 4.0.4 runs but PySAT exposes no statistics for it; CryptoMiniSat's backend is not installed.

165,725 conflicts to *verify* a fully determined 7-stage computation over an 8-bit prime.

### 4c. CNF size scaling, and the extrapolation to the real instance (`runs/cnfsize.json`)

| m | stages | vars | clauses | clauses / stage |
|---|---|---|---|---|
| 8 | 8 | 57,299 | 388,365 | 55,481 |
| 12 | 12 | 191,699 | 1,360,012 | 123,637 |

Both rows fit `clauses/stage ≈ 863 · m²` (867 and 859) and `vars/stage ≈ 111 · m²` (112 and 111)
— the expected quadratic cost of bit-blasted modular multiplication. Extrapolating to the actual
instance (m = 256, 255 stages):

> **≈ 1.4 × 10¹⁰ clauses over ≈ 1.9 × 10⁹ variables.**

In DIMACS that is roughly **600 GB of text**; this machine has 29 GB of disk. The CNF for the real
instance **cannot be written**, let alone solved. This is the concrete form of "beyond current
solvers": not a guess about solver strength, a file-size measurement.

### 4d. CP — OR-Tools CP-SAT (`runs/cpsat.json`), 300 s, 1 worker

| m | selectors pinned | selectors free |
|---|---|---|
| 8 | OPTIMAL 0.2 s | OPTIMAL 0.7 s, 79 conflicts, 379 branches |
| 10 | OPTIMAL 0.0 s | OPTIMAL 4.6 s, 578 conflicts, 2,644 branches |
| 12 | OPTIMAL 0.0 s | **UNKNOWN 300 s**, 19,334 conflicts, 74,176 branches |
| 14 | OPTIMAL 0.0 s | **UNKNOWN 300 s**, 10,741 conflicts, 173,239 branches |
| 16 | OPTIMAL 0.1 s | **UNKNOWN 300 s**, 13,779 conflicts, 192,009 branches |

**CP-SAT is by far the best of the three families** — its integer propagation handles modular
arithmetic natively instead of bit-blasting, and it verifies a pinned instance instantly at every
size tried. But on the free problem it dies between 10 and 12 bits: 0.7 s -> 4.6 s -> fail, i.e.
roughly **×2.6 per bit of prime**, which is *worse than brute force*. At m=12 exhaustive
enumeration of the 2¹¹ configurations takes about 10 ms; CP-SAT had not finished after 300 s —
**about 10⁵ times slower than enumeration on the same instance.**

### 4e. Option 2 from the brief — uninterpreted stage values (`opt2_uf.py`, `runs/opt2.log`)
Declare an uninterpreted sort with a binary operation, give the law only as axioms
(commutativity + associativity), bit-blast only the selectors:

| selectors | axioms | result | time |
|---|---|---|---|
| 8 | assoc + comm | **sat** | 0.01 s |
| 16 | assoc + comm | **sat** | 0.03 s |
| 8 | + all leaves and target distinct | unknown | 64 s |

**The idea is vacuous, and now measured to be.** Without field semantics the axioms admit models
where the operation collapses, so the solver returns SAT instantly with a model that says nothing
about the instance. Adding enough structure to rule that out (distinctness) immediately pushes the
formula into a quantified fragment z3 cannot decide, and the only way to make the answer *mean*
something is to reintroduce the field arithmetic — i.e. option 1. Searching configurations without
field semantics cannot work, because the target is a field element and nothing else constrains it.

### 4f. Option 3 from the brief — CP/MIP over channel choices
CP-SAT is exactly that experiment and its numbers are in 4d. The additional hard wall: CP-SAT's
integer domains are int64, so the real instance's field elements do not fit at all; a limb
decomposition would be needed, and limb-wise multiplication with carries *is* bit-blasting, which
returns us to 4c's 10¹⁰ clauses. No MIP formulation was attempted beyond this, because a
big-M relaxation over a 256-bit modulus has no usable LP relaxation — the feasible set is
2²⁵⁶ isolated points and the continuous relaxation is the whole box.

## 5. WHY the encodings cannot win, stated as a cost, not an opinion
The reduced problem is scalar recovery in a cyclic group of **256-bit prime** order (section 2).
The best generally applicable methods for that are the meet-in-the-middle / cycle-finding family at
**≈ 2^128 group operations**, and prime order rules out any decomposition into smaller subproblems.
No complete search procedure — CDCL, CP, MIP or otherwise — can be faster than the best method for
the problem it encodes. The measurements above show these tools are in fact *slower than
enumeration* at 12 bits, so they are ~2^120 away, not close.

**This is a cost measurement, not a proof of infeasibility, and nothing here claims the instance is
unsatisfiable.** A solution exists: the unique `k` with `k·L_0 = T`.

## 6. Searches actually run on the real 256-bit instance (not siblings)
- **Exhaustive Hamming weight ≤ 6 — NO SOLUTION** (`search_lw.py`, `runs/lowweight6.log`, 108 s).
  Meet-in-the-middle over all 2,796,417 subsets of size ≤ 3 on each side, so every `k` with at most
  6 one-bits was tested. **k has Hamming weight ≥ 7.** (Knob set: all 256 selector bits, no other
  variable; this is a statement about `k`, not about any wire.)
- **BSGS for a small scalar** (`bsgs.py`), 2²² baby steps, covering `k < 2⁴⁴` — running at the time
  of writing, resumable by re-invoking the script; partial and no hit so far.

## 7. Configuration-dependence of the score, measured (`cfgscan.py`, `defect.py`)
Using agent F's chain repair `gs2.solve` (imported read-only) to build *complete* assignments for
chosen selector configurations, then scoring with F's exact evaluator:

| configuration | score | nonzero atoms | equations touched | failing | satisfied inside footprint |
|---|---|---|---|---|---|
| deliverable ON-set {24601, 2081} | 39,005 | 6 | 28 | 28 | 0 |
| any single bit (27 tested, all identical) | 39,013 | 3 | 20 | 20 | 0 |
| three bits {24601, 2081, 47} | 38,985 | 9 | 48 | 48 | 0 |

Two measured facts:
1. **The defect grows linearly in the number of live leaves** — 3 atoms / 20 equations per live
   leaf, and every touched equation fails. Turning on *fewer* leaves is strictly better under this
   repair, and the all-single-bit result is uniform across every bit tested (so it is a property of
   the shape, not of which leaf).
2. **The deliverable's 39,026 does not come from its configuration; it comes from its defect
   placement.** The same configuration repaired by `gs2` scores 39,005. The deliverable puts its
   7 nonzero atoms into 12 equations of which 5 still cancel; `gs2` puts 6 atoms into 28 equations
   of which 0 cancel. **Placement, not configuration, is worth ~21 equations here.**

Knob set for those rows: `gs2.solve`'s knob set = every free input reaching a nonzero residual
atom, minus the frozen set {the four pinned root/target wires, the selected selector bits, the
tree partner flags}. Selector configuration: exactly the bits listed in the table, all others 0.

## 8. What I did NOT do / open
- Did not run the deliverable's own placement optimiser (`s10/lattice3.py` machinery) on a
  single-bit configuration. That is the one experiment that could still move the score: the
  single-bit configurations have a *smaller* defect (3 atoms vs 7) and the notebook's accounting
  rule is about how many of a footprint's equations can be made to cancel. Nobody has priced the
  single-bit footprints with that method.
- Did not finish BSGS; did not attempt weight 7–10 (needs ≈1.7×10⁸ stored points, ~1.4 GB, feasible).

## 9. Pricing the alternative configurations — they lose (`price.py`, `runs/price.json`)
For a configuration's defect footprint, build the (equations touched) × (nonzero atoms) integer
coefficient matrix and compute the maximum number of those equations a **nonzero** atom vector can
annihilate. This ignores integer/mod-p realizability, so it is an *optimistic ceiling*.

| configuration (placement reached by `gs2.solve`) | atoms | eqs touched | max killable | floor failing | ceiling score |
|---|---|---|---|---|---|
| single bit {24601} | 3 | 20 | 7 | 13 | **≤ 39,020** |
| pair {24601, 2081} | 6 | 28 | 17 | 11 | **≤ 39,022** |
| deliverable's own placement (from `NOTEBOOK.md` §Session 10) | 7 | 12 | 5 | 7 | **39,026** |

**SCOPE — this bounds `gs2`'s SUPPORT, not the configuration and not the instance.** `price.py`
fixes the nonzero-atom support to whatever `gs2.solve` lands in and optimises only the *values* on
that support. Which atoms are nonzero is itself a free choice of the repair — the deliverable
deliberately violates GATE atoms to move its defect into a cheap 12-equation footprint. So the
rows above say: *at `gs2`'s support*, single-bit and pair configurations cannot beat 39,026. They
say nothing about supports a placement optimiser could reach. Quoting these numbers as a property
of the configuration, or of the instance, would be wrong.
Knob set: `gs2.solve`'s knob set (every free input reaching a nonzero residual atom) minus the
frozen set {x22162, x30213, x24468, x18956, the selected selector bits, the tree partner flags};
selector configuration: exactly the bits named in each row, all other selectors 0. This prices
*the placement `gs2` lands in*, not every placement that configuration admits — an optimiser that
relocates the defect could do better, and that is the open experiment in section 8.

## 10. CP-SAT's representation ceiling, measured
`opt3_cpsat.py` at m = 31 returns MODEL_INVALID at once: the product domain `p²` exceeds what
CP-SAT will accept. So CP-SAT cannot even *state* the real instance's arithmetic — the wall is at
about 31 bits of modulus, versus the 256 required.

## 11. Verdict for this angle
Automated reasoning against the reduced problem is **measured to be out of reach**, by a margin of
roughly 2^120, and the measurement is not a solver opinion: the CNF for the real instance is
~1.4×10^10 clauses (~600 GB, versus 29 GB of disk), CP-SAT cannot represent 256-bit field
arithmetic at all, and where the tools *can* be run they are ~10^5 times slower than exhaustive
enumeration of the same instance. Agent C's original abandonment of SAT/SMT was right, and it
remains right after the decode — but now for a stated, reproducible reason rather than an
instance-size heuristic.

## 12. Coordinator's contradiction, resolved: the ceiling does NOT bind, and the experiment was run

**Resolution.** `price.py` fixes the nonzero-atom **support** to whatever `gs2.solve` lands in and
optimises only the *values* on that support. Support is itself a free choice of the repair — the
deliverable deliberately breaks upstream atoms to relocate its defect. So §9's ceilings bound
**my repair's support**, not the configuration and not the instance. §8's experiment was live.
Both documents are now re-quoted that way (§9 SCOPE block, `RESUME_R.md` §2/§3).

## 13. TOOLING BUG, worth flagging to the fleet (`crosscheck.py`, `runs/crosscheck.json`)
Agent F's fast evaluator `E` and the lab's `checker.py` **disagree on the deliverable**:

| evaluator | failing on `best/new_instance_partial_39026.json` |
|---|---|
| `checker.py` (ground truth) | **7** — 12231, 12270, 12350, 14584, 18673, 22044, 29125 |
| agent F's `E` | **13** — the same 7 **plus 2554, 6816, 8124, 8680, 9123, 9421** |

`E` is a strict over-approximation there (`chk-only` is empty), so it scores the deliverable
**39,020 instead of 39,026**. The over-report is *assignment-dependent*, not a constant offset:
re-scoring every assignment I generated with `checker.py` (`rescore.py`, `runs/rescore.json`):

| assignment | checker score | E score | E over-reports by |
|---|---|---|---|
| deliverable | **39,026** | 39,020 | **6** |
| gs2 {24601} | 39,013 | 39,013 | 0 |
| gs2 {2081} | 39,013 | 39,013 | 0 |
| gs2 {47} | 39,013 | 39,013 | 0 |
| gs2 {24601, 2081} | 39,005 | 39,005 | 0 |

So my configuration scores were right, and the "placement is worth 21 equations" figure
(39,026 vs 39,005 on the same configuration) is right **on the checker's scale**. But anyone
scoring with `E` instead of `checker.py` will under-report the deliverable by 6, and the
deliverable's own defect footprint cannot be read off `E` at all.

Corrected footprint of the deliverable: its defect occupies **13 equations of which 6 cancel**
(7 fail). Not the "7 atoms / 12 equations / 5 cancelling" I quoted earlier — that came from
`NOTEBOOK.md` §Session 10, whose atom numbering (22229, 35758…) is a *different indexing* from
`E`'s (3130, 7251…). Those two numberings must not be cross-quoted; I did, and it was wrong.

## 14. THE EXPERIMENT — why a single-bit configuration cannot be given cancellation (`cancel.py`)
The instrument is cancellation, not support. Measured for the single-bit footprint
(3 live atoms, 20 equations, **0 cancelling**), against the deliverable's 13 equations with
**6 cancelling**:

| | count |
|---|---|
| equations in the single-bit footprint | 20 |
| ...that can *never* cancel (only one atom in the whole equation) | **0** |
| ...that do have a dead partner atom available | **20** |
| ...whose cheapest available partner touches ≤ 3 other equations | **0** |

So cancellation is structurally *available* everywhere, and expensive everywhere.
**The cheapest partner atom anywhere in the footprint is atom 7954, which occurs in 10 other
equations**; next are 4490 and 4497 at 11, then 4496 and 4561 at 12, 8259 at 13, 4500 and 4331
at 14. The single-bit configuration is 13 equations short of beating the deliverable
(39,013 against a 39,027 target), so it needs ~7 purchased cancellations.

**Stated at the strength the measurement supports:** turning on a partner to cancel one equation
**touches at least 10 further equations**. *Touched is not failed* — those 10 may themselves
cancel downstream, and this lab has three independent demonstrations that incidence does not price
cost (agent C check-in 3; agent P's withdrawn carrier check-in 10; agent L measuring the incidence
scorer inflating the deliverable's own cost). So this is **not** a proof that the purchase cannot
pay; it is a measured price floor in *touches*, with a ~10x margin against the 13-equation deficit.
The margin is what makes the conclusion likely, not the arithmetic.

Scope: the partner-occurrence counts are facts about `EQUATIONS.txt`'s incidence structure and are
not knob-set dependent. The *inference from touches to failures* is not established, and which
footprint `gs2` lands in remains scoped to my repair.

**Why the deliverable wins, stated cleanly:** it is not using a better configuration and not a
larger support. It sits in a rare footprint where 6 of 13 equations cancel *for free* — no partner
atoms had to be bought. Every footprint I reached charges ≥10 equations for the first cancellation.

## 15. INVERTED SEARCH — rank footprints first, ask reachability second
(`footprints.py`, `reach2.py`, `rank.py`, `relax.py`; raw in `runs/footprints1.json`,
`runs/reach2.json`, `runs/rank.json`, `runs/relax.json`.)

Every search in this campaign went configuration-first. This goes the other way: for a live-atom
support S define
`cost(S) = |equations touched by S| − max #{touched equations a nonzero value vector annihilates}`,
which is the failing floor for **any** assignment whose nonzero-atom set is exactly S. Rank by
cost, then ask which configurations route into the cheap ones. Beating 39,026 needs cost ≤ 6.

E's model: **9,032 atoms, 39,033 equations, 121,261 incidences.**

### 15.1 The cost distribution, and why the deliverable is already at an optimum
Single-atom footprint costs (cost -> how many atoms): 2->1, 3->1, 4->1, 5->1, 6->2, 7->8, 8->36,
9->113, 10->354, 11->863, 12->1481, 13->1793, 14->1835, 15->1327, 16->721, 17->339, 18->112,
19->33, 20->9, 21->1.

Splitting those by kind is the whole result:

| | min single-atom cost |
|---|---|
| **relational** atoms (≥2 variables — the ones that can carry the defect) | **7** |
| **boolean-ness** atoms `x·(1−x)` (one variable, assert `x ∈ {0,1}`) | **2** |

**Every atom with cost ≤ 6 is a boolean-ness atom on a single variable.** The cheapest relational
atoms cost exactly 7 — 3131, 3138, 3588, 8749, 8777 — and **atom 3131 is one of the deliverable's
own four live atoms.** So the deliverable sits on the cheapest relational footprint that exists,
and **7 failing is the single-atom optimum over relational supports.** That is why five agents'
configuration-first searches all bottomed out at the same number.

`|S| = 2` does not help either (`rank.py` step 3): over the 60 cheapest relational atoms, all 34
equation-sharing pairs floor at **≥ 7**, the best (3131+3138, 3138+8749, 3588+5558, 8777+2569)
tying at exactly 7 (9 touched, 2 killed). Cancellation buys back exactly what the larger support
costs.

### 15.2 The cheapest atom of all is a trap, and the trap is instructive
Atom 8508 = `x29570·(1−x29570)`, cost **2**. But `reach2.py`: x29570 occurs in **no other atom**.
It is cheap *because it is disconnected*, and a disconnected variable cannot carry the defect —
the defect is the mismatch between the leaf selection and the target and must flow through atoms
that link them. Same for 7888 (`x27026`, cost 5) and 8512, 8764 — all zero-collateral.
**Cheapness and load-bearingness are anti-correlated here.** This is the ranking's own version of
the incidence-does-not-price-cost trap, and it eliminates the four cheapest entries.

### 15.3 THE LEVER: relax a SELECTOR off {0,1}
Of E's 2,283 boolean-ness atoms, **173 are on selector / conditional-pin variables**, and the
cheapest of those are *not* disconnected:

| atom | variable | cost | selector? |
|---|---|---|---|
| 8509 | x33095 | **3** | yes (tree A) |
| 7889 | x19326 | **6** | yes (tree A) |
| 8510 | x28825 | **6** | yes (tree A) |
| 8511 | x4362 | 7 | yes (tree A) |

A selector `b` off `{0,1}` does **not** force its mux atoms nonzero: the mux relation
`acc' = acc + b·(S − acc)` is still satisfiable — `acc'` just becomes a free point on the line
through `acc` and `S`. So the collateral can stay zero and the only atoms forced nonzero are the
boolean-ness atoms of the relaxed selectors. (`reach2.py`'s "true-touch 18–22" for these assumed
all collateral goes nonzero; it does not have to, which is why that number is an over-estimate.)

Counting parameters: one relaxed selector is 1 free field parameter against the 2 coordinates of
the target — generically no solution. **Two relaxed selectors are 2 parameters against 2
conditions, generically solvable for any boolean choice of the other 254.** Pricing the pairs by
the union of their boolean-atom equations:

| relaxed pair | equations in union | overlap | floor |
|---|---|---|---|
| **x33095 + x19326** | **6** | 3 | **39,027** |
| **x33095 + x28825** | **6** | 3 | **39,027** |
| x33095 + x4362 | 7 | 3 | 39,026 (ties) |
| x19326 + x28825 | 7 | 5 | 39,026 (ties) |
| x28825 + x4362 | 7 | 6 | 39,026 (ties) |
| x19326 + x4362 | 8 | 5 | 39,025 |

**Two pairs floor at 6 failing equations — 39,027.** The overlap is what does it: x33095's 3
equations are a subset-heavy intersection with x19326's 6, so the union is 6 rather than 9.

### 15.4 Status of the 39,027 target — UNREALIZED, and it is a floor, not a score
`cost(S)` is a lower bound on failing equations for a support; it says nothing about whether an
assignment with that support exists. Realizing it needs the fold to reach the target exactly with
two non-boolean coefficients — a 2-unknown solve whose difficulty I have not measured, and the
composite map is a rational function of both parameters through every downstream stage, so its
degree may be enormous. `realize.py` hands the two relaxed selectors to `gs2`'s repair at fixed
non-boolean values and scores with `checker.py`; that is a weak attempt (gs2 propagates forward
from the leaves, so it cannot back-solve the two parameters against the root). **It failed:** with
either relaxed selector frozen at a non-boolean value `gs2.solve` returned nothing in ~15 minutes
per call (it repairs forward and appears to loop trying to restore booleanness), so not one
candidate was scored. Stopped, cores released.
**No score above 39,026 has been verified. 39,027 is a floor with an unrealized construction, not
a result, and must not be quoted as a score.**

What realizing it actually needs: a *backward* solve, not a repair. Fix the 254 boolean selectors,
treat the two relaxed ones as unknowns `t1, t2`, push the accumulator recurrence
`acc' = acc + t*(S - acc)` symbolically to the root, and solve the resulting 2x2 system against the
target's two coordinates. The obstruction to price FIRST: every stage downstream of a relaxed
selector applies the chord law to an off-curve point, so the composite degree grows with the number
of stages after it. **If x33095 / x19326 / x28825 sit near the root the system is small and this is
straightforward; near the leaves it is hopeless. Nobody has checked which. That is the single next
measurement and it is one cheap lookup.**

## 16. Depth profile, an obstruction I priced — and then falsified myself

### 16.1 The lookup, fixed and validated (`depth.py`)
`agentF_work/pins.json` stores each pin as `[[wire, val], [wire, val]]` where the **second** entry
is the x-coordinate in the unshifted frame — the pair is (y, x), not (x, y). So
`ladder index = LX[(val_second + S) mod P]`. Validated against the four-agent-confirmed ON-set:
**x24601 -> 72 OK, x2081 -> 235 OK.**

Positions of the cheap relaxable selectors:

| atom | var | cost | ladder index | stages downstream |
|---|---|---|---|---|
| 7887 | x24267 | 4 | **8** | 247 |
| 8509 | x33095 | 3 | **132** | 123 |
| 7889 | x19326 | 6 | **73** | 182 |
| 8510 | x28825 | 6 | **218** | 37 |
| 8511 | x4362 | 7 | **243** | 12 |

### 16.2 A better floor my earlier scan missed
`relax.py` only ranked the 25 cheapest selector-boolean atoms and a regex bug dropped atom 7887.
The exhaustive scan over all **253** placed selector-boolean atoms (`obstruct.py`) finds

> **x24267 (ladder 8) + x33095 (ladder 132): union of only 4 equations -> floor 39,029.**

Four pairs beat 39,026 (floors 39,029, 39,027 x3); four more tie.

### 16.3 The obstruction I computed (`tradeoff.py`) …
What drives the backward solve is the **gap** between the two relaxed ladder indices, not depth
from the root: stages after the deeper selector are boolean and invertible (walk the root back to
a concrete required point, no degree growth), but each stage *between* the two relaxed selectors
applies the chord law to an off-curve accumulator and roughly doubles the degree in `t1`, so the
elimination is degree ~2^gap. Best floor reachable at each gap budget:

| gap ≤ | best floor |
|---|---|
| 1, 2 | 39,019 |
| 4–16 | 39,021 |
| 24 | 39,024 |
| 32–48 | 39,026 (ties, never beats) |
| **59** | **39,027 — first beat** |
| 128 | **39,029** |

Monotone, no cheap corner: every pair that beats 39,026 needs gap ≥ **59**, i.e. elimination
degree ≥ 2^59. On that reading the lever is priced out.

### 16.4 … and why that obstruction is WRONG (`collapse.py`)
**The 2^gap figure assumes the selectors between the two relaxed ones are arbitrary. They are not
— I choose them.** Set every selector strictly between (and after) the relaxed pair to **0** and
the mux `acc' = acc + b·(S − acc)` becomes the **identity**: the accumulator does not move, no
chord is applied to it, and the degree does not grow at all. The gap becomes irrelevant.

Tested on scaled siblings of identical shape, brute-forcing `t1` over the small field
(`collapse.py`, m = 8…16, gaps 1…14): **solutions exist independently of gap** — 10 of 15 cases
had a `(t1, t2)` reaching the target, with gap 1 sometimes failing and gap 14 succeeding. The hit
rate is what a low-degree elimination with O(1) expected roots predicts, not a 2^gap wall.

**I priced an obstruction and then falsified it with my own experiment. §16.3's table is
superseded; it is retained because the reasoning error in it is the instructive part.**

### 16.5 The solve, at the REAL 256-bit prime (`solve2.py`)
With base accumulator `A` and relaxed leaves `i < j`, everything else off:
`acc = A + t1·(chord(A,L_i) − A)` (degree 1), and `T` must lie on the line through `acc` and
`chord(acc,L_j)`. Interpolating that residual and root-finding by `gcd(t^P − t, f)`:

| relaxed leaves | floor | elimination degree | deg gcd(t^P−t, f) | |
|---|---|---|---|---|
| **8 + 132** | **39,029** | **8** | 1 | **roots exist** |
| 73 + 132 | 39,027 | 8 | 3 | roots exist |
| 8 + 73 | 39,027 | 8 | 1 | roots exist |
| 132 + 218 | 39,027 | 8 | 2 | roots exist |

**Degree 8, not 2^124.** The two-parameter system is solvable at the real prime for all four
beating pairs.

### 16.6 Status — still NOT a score
What exists: a floor of 39,029, and a field solution `(t1, t2)` for the parameters. What does
**not** exist: a materialised 38,748-wire assignment, and any `checker.py` verification.
Two things must be checked before this is worth anything:
1. `solve2.py` uses my **ladder-chain model with the accumulator seeded at `L_0`**. The real
   circuit is a *tree*; if its accumulator base or gating differs, `A` is wrong and the roots do
   not transfer. **This is the load-bearing assumption and it is unverified.**
2. Even with correct `(t1, t2)`, every other wire must be filled consistently. `gs2` cannot do it
   (it repairs forward and restores booleanness); a forward evaluator that accepts non-boolean
   selectors is needed.

**Nothing above 39,026 is verified. 39,029 is a floor plus a parameter solve, not a score.**

## 17. STEP 1 — validating `A` against the deliverable. IT FAILS. I am withdrawing the roots.

### 17.1 What the deliverable actually contains (`validate_A.py`, `validate_A2.py`, `validate_A3.py`)
It specifies **3,540 of 38,748** variables, **94 distinct values** after reduction mod P.

| check | result |
|---|---|
| pin/selector variables nonzero in the deliverable | **exactly 2 of 256: x2081 and x24601** — my ON-set |
| the 4 coordinate wires those 2 pins name | **all 4 set, all 4 holding exactly the value `pins.json` names** |
| the other 505 coordinate wires | unset |
| ladder points with a coordinate on a wire (reduced mod P) | **2 of 256 — leaves 72 and 235, i.e. exactly the live ones** |
| target coordinates on wires | **yes** (x on 13682/22162/24468, y on 10156/18956/30213) |
| **leaf 0** (my accumulator seed) | **absent** |
| **fold `L72 + L235`**, any frame | **absent** |
| `L0 + L72`, `L0 + L72 + L235` | **absent** |

(My first pass reported 0/256 ladder points on wires. That was **my own bug** — wire values are
unreduced ~89-digit integers and I searched for the reduced/shifted coordinate. Same class of
mistake as the pins lookup. Corrected above; the corrected numbers are the ones that matter.)

### 17.2 The verdict
Two separate things came out of this, and they point opposite ways:

**Confirmed:** leaf gating is real and my ON-set reading is right. Pin = 1 puts that leaf's two
coordinates on the two wires `pins.json` names, exactly, and only the two live pins do so. That is
a direct check of x24601 -> 72 and x2081 -> 235 against wire contents, stronger than my earlier
`fold(k) != T` test — which, I should say plainly, was nearly **vacuous**: almost any wrong model
also yields `fold != T`.

**Refuted:** the accumulator. `solve2.py` seeds the chain at `L_0`. **Leaf 0 is absent from the
deliverable while both live leaves are present** — so if the fold were seeded at `L_0`, the
deliverable's live-leaf set would be {0, 72, 235} and the ON-set would carry a `2^0` term. It does
not; four agents independently read it as `2^72 + 2^235`. **`A = L_0` is wrong.**

And no accumulator value of any kind appears: not the fold of the two live leaves, not any partial.
The deliverable holds *the inputs and the target and nothing between them* — consistent with T's
and Q's finding that routing is a constraint that is not propagated.

That leaves my model unable to be checked against the only verified object in the lab: it predicts
wire contents that the deliverable simply does not contain. **A model that makes no prediction the
one verified point can test cannot be trusted to produce a new one.**

### 17.3 What I am withdrawing, and what survives
**WITHDRAWN — do not carry forward:**
* the four `(t1, t2)` root results of §16.5, including the 39,029 pair. They are roots of a system
  built on `A = L_0`, and that seed is refuted. **Do not route anything to agent M.** Step 2 is
  not reached, by my own step-1 test.
* the degree-collapse argument of §16.4 *as an instance claim*. It is still true of the siblings,
  but the siblings are my own construction and the collapse depends on the same chain model.
* the claim that relaxing a selector leaves the mux atoms satisfiable — that rests on the mux form
  `acc' = acc + b·(S − acc)`, which is model, not measurement.

**SURVIVES — pure equation-incidence, no group model involved:**
* §15.1: every atom with footprint cost ≤ 6 is a boolean-ness atom; **minimum cost over relational
  atoms is 7**; atom 3131 is one of the deliverable's own live atoms; `|S| = 2` ties at 7 but never
  beats. This is the explanation for why configuration-first searches bottom out at 7 and it does
  not depend on any of the above.
* §15.2: the cheapest atom (8508) is cheap because disconnected — cheapness and load-bearingness
  anti-correlate.
* the floor **39,029** for `{x24267, x33095}` as a *4-equation union* — an incidence fact. What
  dies is the argument that it is reachable, not the arithmetic that it is a floor.
* §13: the `E` vs `checker.py` discrepancy, and §16.1's corrected pins lookup (validated against
  x24601 -> 72, x2081 -> 235 by wire contents).

### 17.4 One measurement that disagrees with a relayed claim — reporting, not adjudicating
I was told T and Q measured that "the deliverable keeps its live leaves with all 256 selectors
forced to 0". **I measure exactly 2 of the 256 `pins.json` variables nonzero in the deliverable
(x2081, x24601 = 1), and their 4 named coordinate wires all hold the named values.** Either their
"selector" is a different variable family than F's `pins.json`, or one of us is wrong. I have not
read their directories and am not adjudicating — flagging it so someone who can, does. My numbers
are reproducible from `validate_A2.py` in three lines.
