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
