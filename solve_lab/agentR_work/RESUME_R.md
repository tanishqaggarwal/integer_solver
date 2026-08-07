# RESUME_R — agent R (automated reasoning against the REDUCED problem). Self-contained.

## 0. Score
- Baseline re-verified by me from cold, 26 s:
  `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
- **I did not beat it.** No assignment of mine scores above 39,013.
- **No infeasibility claim is made anywhere in this directory.** A solution exists.

## 1. MAIN RESULT — the reduced problem, made explicit
Re-derived here from agent F's decode artifacts (read-only), every step checked:
`model.py` `group.py` `ladder.py` `order.py`.

1. `X = x + K/3` removes the offset: the stage law becomes the plain chord law (checked vs
   `chordK` on 200 random pairs).
2. All 248 forced pin points **and the target** satisfy one relation `Y^2 = X^3 + B`,
   `B = 64019533680030876408443198762210829058751700634554282185987325820393598524794`
   (fitted from 2, verified on 246 + target, 0 exceptions).
3. So the fold is an **abelian group law** — commutativity and associativity checked on 200 random
   triples. **The fold of a leaf subset therefore does not depend on the tree shape.** F's
   outstanding decode (56 slot pairs, 24 leaf-adjacent stages) is *not needed*.
4. **The 256 leaves are one doubling chain** `L_i = 2^i · L_0`: 9 doubling-closed pieces
   (111,61,28,12,11,9,8,5,3 = 248 named), each piece's end doubled twice lands exactly on the next
   piece's start (8 splices, 1 unknown point per gap, 1 head) -> `ladder.json`, 256 points.
5. **Group order** by Cornacchia (`4p = L²+27M²`), the candidate annihilating the base point:
   `N = 115792089237316195423570985008687907852837564279074904382605163141518161494337`,
   **256-bit prime**. Prime order ⇒ no decomposition into smaller subproblems.

> **The reduced problem is: find the 256-bit scalar `k` with `k·L_0 = T`,
> bit *i* of `k` being the selector of leaf `L_i`.**

**F's requested validation PASSED**: deliverable ON-set {24601, 2081} -> ladder indices {72, 235}
-> `k = 2^72 + 2^235`, and `fold(k) != T`, exactly as F said a correct evaluator must report.

## 2. Measured outcomes (full tables in `LOG.md` §4, raw in `runs/`)
- **z3 QF_BV**, m=8 sibling (7 stages, 8-bit prime), *all selectors pinned to the known solution*:
  **sat in 119 s**. Selectors free: **timeout 300 s**. QF_NIA: timeout even pinned.
- **CNF** (z3 bit-blast, m=8 pinned, 57,299 vars / 388,439 clauses):
  CaDiCaL 1.9.5 sat 107.6 s, **165,725 conflicts**, 354,424 decisions, 34.6M propagations;
  Glucose 4.2 sat 276.7 s, 353,100 conflicts, 122.9M propagations.
- **CNF size scaling**: clauses/stage = 55,481 (m=8) and 123,637 (m=12), both fitting
  `≈863·m²`. Extrapolated to the real instance (m=256, 255 stages):
  **≈1.4×10^10 clauses / 1.9×10^9 vars ≈ 600 GB of DIMACS**, versus 29 GB of disk.
- **CP-SAT** (best of the three families): free selectors solved at m=8 (0.7 s, 79 conflicts) and
  m=10 (4.6 s, 578 conflicts), then **UNKNOWN at 300 s for m = 12, 14, 16, 20, 24, 28**.
  At m=12 brute force takes ~10 ms, so CP-SAT is ~10^5× slower than enumeration.
  At m=31 CP-SAT returns MODEL_INVALID — it cannot represent the arithmetic at all.
- **Uninterpreted-law variant** (brief's option 2): sat in 0.01 s with no field semantics —
  **vacuous**; adding distinctness makes it `unknown` in a quantified fragment. Closed.
- **Exhaustive Hamming weight ≤ 6 on the real 256-bit instance: NO SOLUTION** (108 s,
  meet-in-the-middle over all size-≤3 subsets both sides). So `k` has weight ≥ 7.
- **Configuration dependence**: single-bit configurations give 3 nonzero atoms / 20 equations
  (uniform over 45 bits tested), pairs 6/28, triples 9/48; the deliverable's *placement*
  (7 atoms / 12 equations, 5 cancelling) is worth ~21 equations more than `gs2`'s repair of the
  same configuration.
  **SCOPED CEILING — read the scope before quoting the number.** `price.py` fixes the nonzero-atom
  *support* to the one `gs2.solve` happens to land in and then optimises only over the *values* on
  that support. Under that fixed support the ceilings are single bit **≤ 39,020** and pair
  **≤ 39,022**. **This is a statement about `gs2`'s support, not about the configuration and not
  about the instance.** A repair that relocates the defect to a different support is not bounded by
  it — and relocating the defect is exactly what the deliverable does. The ceiling therefore does
  **not** kill the §4 experiment.

## 3. Confirmed / refuted
- CONFIRMED: F's tree decode, its law, its invertibility, and its ON-set prediction.
- CONFIRMED (new, unconditional): the fold is a commutative associative group law, so the tree
  shape is irrelevant and the leaves form a single doubling ladder over a 256-bit prime order.
- REFUTED: that the reduced problem is 96 coupled stage constraints — it is one scalar equation.
- REFUTED: that "encode the reduced problem in SAT/SMT/CP" is a live attack. Measured, not assumed.
- NOT REFUTED, and previously mis-stated by me: whether a different selector configuration can
  beat 39,026. What is measured is only that **`gs2.solve`'s support** for single-bit and pair
  configurations is capped at 39,020 / 39,022. Support is a free choice of the repair, so this
  bounds my repair, not the configuration and not the instance. See §4.

## 4. Highest-EV next experiment
Run the deliverable's own placement optimiser (`solve_lab/s10/lattice3.py` machinery: integer
kernel by column HNF + the two-congruence integer system) on a **single-bit** configuration. Its
defect is smaller (3 atoms vs 7) but `gs2` places it badly (20 equations, 0 cancelling). Nobody has
priced the single-bit footprints with the method that produced 39,026. If a single-bit footprint
admits a placement into ≤ 12 equations with ≥ 6 cancelling, that beats the deliverable.
Second: extend the weight search to 7–8 (~1.7×10^8 stored points, ~1.4 GB — feasible);
`bsgs.py` (k < 2^44) was still running at handoff and is resumable by re-invoking it.

## 5. Files
Model/derivation: `model.py group.py ladder.py order.py fastgrp.py ladder.json points_short.json`
Encoders: `z3enc.py encode.py witness.py sibling.py` -> `encodings/*.smt2`, `encodings/*.cnf`
Benchmarks: `bench.py cnfbench.py opt3_cpsat.py opt2_uf.py` -> `runs/bench.json runs/cnf.json
runs/cnfsize.json runs/cpsat.json runs/opt2.log`
Real-instance searches: `search_lw.py bsgs.py cfgscan.py defect.py price.py` ->
`runs/lowweight6.json runs/bsgs.json runs/cfgscan.json runs/defect.json runs/price.json`
Narrative with every number: `LOG.md`.

## 6. Standing rules honoured
- Every "nothing can move X"-style statement above carries its knob set and selector configuration
  (LOG.md §7, §9).
- No state I produced exceeds 4,300 digits, so `checker.py` was valid for everything I verified;
  `agentE_work/verifyE.py` was not needed.
- No generator forensics, no named-curve framing: `B`, the ladder and `N` are all *measured* from
  the decoded law and reproducible by the four scripts named in §1.
