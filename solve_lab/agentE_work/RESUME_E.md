# RESUME_E — agent E, self-contained handoff

## 0. Two rules that apply to everything below
**(a) Verification.** My states contain integers up to ~4,430 decimal digits, above Python's
4,300-digit string cap, so **`checker.py` cannot PARSE them** (ValueError, not a wrong answer).
Use `python3 verifyE.py <file>` — it raises only that cap and then calls
`checker.load_equations`, `checker.load_assignment`, `checker.evaluate_all` **unmodified**.
`checker.py` itself is untouched.  Always say this explicitly; "checker.py says" would be false.
**(b) Standing caveat (earned the hard way — four retractions, LOG 14/15, 20/21, 23).**
Every mod-p rigidity / pinning / divisibility statement in this lab is conditional on
(i) the **selector configuration** and (ii) the **knob set** it was measured over.  Both the
coefficients AND the targets move with configuration.  Re-quote any such claim with both, or
it is not a claim about the instance.  My four false barriers were each a property of a
filtered knob set reported as a property of the instance: booleans excluded; affine-only;
a residue multiset assumed fixed; dual-reaching knobs only.

## 1. Scores
- Baseline deliverable `../best/new_instance_partial_39026.json` = 39026/39033 (parses fine
  with plain `checker.py`).  CONFIRMED.
- **My best: 39,005/39,033 — `triple8_39005.json`** (verify with `verifyE.py`).  Only **two**
  nonzero atoms in the whole instance (a20215, a28647) vs the deliverable's eight.
  Seed: `triple8_seed.json` (a dict {var: value} of free-input overrides).

## 2. Rebuild from scratch (~1 min)
    python3 parse3.py       # -> model3.pkl : 39,033 eqs = outer * Z-combo of 40,727 ATOMS
    python3 dag.py          # -> dag.pkl    : 35,004 defs x_out - RHS, ACYCLIC, 8,365 free vars
    python3 -c "import harness"   # -> orient.pkl (orientation + topological order)
    python3 prop2.py        # free=0 propagation -> 38,998/39,033, only 3 violated atoms
Modules: `engine.py` (forward map, cone eval, exact single-var solve), `fast.py` (incremental
downstream-only re-evaluation, 0.08 s/probe, verified exact), `sparse.py` + `intsolve.py`
(singleton/unit-pivot elimination then HNF), `iterfix.py` (iterated closure+solve),
`channels.py` (channel measurement + LOG-16 simultaneous solve), `chanenum.py`, `verifyE.py`.

## 3. How 39,005 was reached (each step exact, re-verified)
1. Seed all free vars 0, propagate -> 3 violated atoms.  Core residual = `OR(a,b)=1` forced
   (a=x_7715, b=x_34554, OR-trees over 178/78 free bits), a 2-way MUX (a20212), and a20215.
2. (1,1) branch closes all four core atoms: one free bit on in each tree, then
   `x_22162 = x_13682`, `x_30213 = x_18956 - x_32237` (4-iteration fixpoint).
3. Remaining = the two bits' pin atoms `b*(free-K) = m*handle`, `pin = p*handle`.
4. The size>=2 obstruction triple `{722,724,726}` **closes exactly** (`triple4.py`): its three
   conditions are really two (p|U plus the exact row forces p|V, gcd(15322661,p)=1), and (U,V)
   is exactly affine in the two NON-BOOLEAN knobs **x_30468, x_33169**.  Closed form, one free
   parameter k, second congruence modulus 1 so every k works; handles x_34496, x_3193 finish it.
5. `triple8.py` then solves 44 of 46 rows of the affine system simultaneously -> **39,005**.

## 4. What binds now — two rows, and the correct model of them
    a20215 : x_24530 - x_5647*x_24908
    a28647 : x_36433 - (x_36990 + x_19239)
Reducing the 5-row cluster {7389,10187,20212,20215,28647} mod p leaves two congruences in the
selector bits.  The reachable directions have **RANK 2** (unit directions x_14853, x_31339),
so it is NOT a rank obstruction.

## 5. SATURATION — the central mechanism (LOG 23)
Per-bit deltas **do not add, they saturate**.  Measured: 2 or 3 bits of one class give exactly
the delta of ONE bit.  A class contributes its coefficient at most once — a channel that is
live or not.  This retroactively explains every additivity failure I recorded (LOG 10, 15.1,
20).  Independently reproduces the pass-through law the fleet derived from circuit structure.
**Consequence: never model multiple selector flips as a sum.**

## 6. The channel model, and the partition (LOG 24-25)
At cfg0 the 256 booleans in the cluster cone split into **3 pairwise-disjoint channels of
178 / 41 / 21**, plus **16 that move nothing**.  `41 + 21 + 16 = 78`, so the cone splits
**178 | 78** — matching the reported root-gate slot supports 178 and 78, reached independently
from the residual side.  Class count and targets are configuration-dependent
(cfg5: 4 channels; cfg7: 4 channels of 41/24/9/4).
`chan_cfg0.json` holds the exact channel membership lists.

## 7. Channel enumeration: RUN, complete at cfg0, negative
All 2^3 channel-sets x 2 representatives = 27 exact (channel-set, representative) pairs through
the simultaneous solve.  **The empty set wins at 39,005**; one live channel 38,969-38,992, two
38,934-38,959, three 38,872-38,913 — monotone in the number of live channels.  cfg5/cfg7
(4 channels) show the same monotonicity.  Command: `python3 chanenum.py <n_representatives>`.

## 8. Next experiments, in priority order
1. **Sweep representatives properly.** Only 2 per channel were tried; the representative
   determines its own pin rows and the 106 bits with known-solvable pin systems
   (`bitsol_*.json`, `scan_B.pkl`, `scanfork_A.pkl`) are the ones to use first.  Cheap: each
   pair is ~3 s.
2. **Test the 16 inert booleans.** They are predicted to be a branch that is not live at cfg0.
   Find a configuration that makes them move; if one exists it is a new channel and the
   partition claim in section 6 becomes checkable rather than numerological.
3. **Monotonicity is suspicious.** Every live channel costs, at every configuration measured.
   If the intended witness has channels live, then either the cost is recovered by pin repairs
   my closure does not reach (raise `maxr`/`maxv` in `channels.simsolve`), or the cluster must
   be attacked at a configuration where a20215/a28647 are not the residual at all.
4. Undecided singletons: 8 with over-budget HNF cores, 6 timeouts, in `runs/scanA6.log`.
   Raise `sparse.solve_sparse(..., maxcore=, maxcorebits=)`.

## 9. Complete scan results (for reuse)
Per-bit pin feasibility, exact: **a-tree 56/178 feasible, b-tree 50/78 feasible**; all 106
verified with zero residual atoms outside the selector core, each giving 39,017
(`bitsol_<bit>_39017.json`).  Infeasible ones fail on an explicit p-divisibility row.
