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

## 4. The single-bit experiment: RUN, and it fails for a quantified reason
Cancellation, not support, is the instrument (`cancel.py`, `runs/cancel.json`). Single-bit
footprint = 3 live atoms, 20 equations, **0 cancelling**; deliverable = 13 equations, **6
cancelling**. Of the 20, **0 can never cancel** (all have a dead partner atom) but **0 have a
partner costing ≤3**: the cheapest partner anywhere is atom 7954, present in **10** other
equations, then 4490/4497 at 11, 4496/4561 at 12, 8259 at 13, 4500/4331 at 14. Each cancellation
bought **touches ≥10 further equations**; the configuration is 13 short.
**Touched is not failed** — those equations may cancel downstream, and incidence has been shown
three times in this lab not to price cost. So this is a measured price floor in *touches* with a
~10x margin, **not** a proof the purchase cannot pay. The occurrence counts are facts about
`EQUATIONS.txt`; the touches-to-failures inference is not established; which footprint `gs2`
reaches stays scoped to my repair.

**Why the deliverable wins:** not a better configuration, not a bigger support — a rare footprint
where 6 of 13 equations cancel *for free*. Every footprint I reached charges ≥10 for the first.

## 4b. INVERTED SEARCH - footprints ranked first, reachability second (LOG.md 15)
E's model: 9,032 atoms / 39,033 equations / 121,261 incidences. Footprint cost
`= |equations touched| - max killable by a nonzero value vector` = the failing floor for a support.

- **Min cost over RELATIONAL atoms (>=2 vars, the only ones that can carry the defect) = 7.**
  Cheapest are 3131, 3138, 3588, 8749, 8777 - and **3131 is one of the deliverable's own live
  atoms.** So the deliverable sits at the single-atom optimum; `|S|=2` ties but never beats it (all
  34 equation-sharing pairs among the 60 cheapest relational atoms floor at >=7).
  **That is why five configuration-first searches all bottomed out at 7.**
- **Every atom with cost <=6 is a boolean-ness atom `x*(1-x)` on one variable.** The cheapest of
  all (8508, cost 2) is a trap: its variable occurs in no other atom - cheap *because
  disconnected*, and a disconnected variable cannot carry the defect.
- **THE LEVER (new): relax a SELECTOR off {0,1}.** 173 of the 2,283 boolean-ness atoms sit on
  selector/conditional-pin variables and the cheap ones are not disconnected: x33095 (cost 3),
  x19326 (6), x28825 (6), x4362 (7). A non-boolean `b` does **not** force the mux atoms nonzero -
  `acc' = acc + b*(S-acc)` stays satisfiable with `acc'` a free point on a line - so only the
  boolean-ness atoms are forced. Two relaxed selectors = 2 free parameters against the target's 2
  coordinates, generically solvable for *any* boolean choice of the other 254.
  **x33095+x19326 and x33095+x28825 each have a 6-equation union (overlap 3) -> floor 39,027.**
- **UNREALIZED.** `realize.py` (gs2 repair, selector frozen non-boolean) returned nothing in
  ~15 min/call - it repairs forward and cannot back-solve two parameters against the root.
  **39,027 is a floor with no construction. Not a score. Do not quote it as one.**

## 4c. THE NEXT MEASUREMENT, and it is cheap
Realizing 39,027 needs a *backward* solve: fix 254 boolean selectors, treat `t1,t2` as unknowns,
push `acc' = acc + t*(S-acc)` symbolically to the root, solve 2x2 against the target. Price the
obstruction FIRST: every stage after a relaxed selector applies the chord law to an off-curve
point, so composite degree grows with the number of downstream stages.
**Find where x33095, x19326, x28825 sit in the ladder. Near the root -> small system, straight-
forward. Near the leaves -> hopeless. Nobody has checked.**
I tried the lookup and it did not converge in the time left: `agentF_work/pins.json` maps a pin
variable to `[[wire_x, val_x], [wire_y, val_y]]`, and feeding `(val_x, val_y)` through
`model.to_short` then indexing `ladder.json` returns NOT FOUND *even for x24601 / x2081*, which I
had previously mapped to ladder indices 72 / 235 by a different route. So the naive reconstruction
above is wrong somewhere, not the ladder. Redo it the way the ON-set validation did (that path
worked); do not trust the two-line version.

Also: extend the weight search to 7–8 (~1.7×10^8 stored points, ~1.4 GB — feasible);
`bsgs.py` (k < 2^44) was still running at handoff and is resumable by re-invoking it.

## 4c. TOOLING BUG — flag to the fleet
Agent F's evaluator `E` reports **13** failing on the deliverable where `checker.py` reports **7**
(the 7 real ones plus 2554, 6816, 8124, 8680, 9123, 9421), scoring it 39,020 instead of 39,026.
The over-report is assignment-dependent, not a constant: on all five `gs2` assignments I generated
`E` and `checker.py` agree exactly. So scores computed with `E` are usable, but **the deliverable's
own footprint cannot be read off `E`**, and anything compared against it via `E` is off by 6.
Also: `NOTEBOOK.md` §Session 10's atom numbering (22229, 35758…) is a different indexing from `E`'s
(3130, 7251…) — I cross-quoted them once and it was wrong; corrected in LOG.md §13.

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
