# Agent A — RESUME (exact integer linear algebra / lattice angle)

## Best verified score: 39,026  — file: solve_lab/agentA_work/A_best_39026.json
(byte-identical to the lab baseline; re-verified with solve_lab/checker.py -> 39026/39033,
failing [12231,12270,12350,14584,18673,22044,29125]).  I did not beat it.

## Re-entry
```
cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py   # caches are NOT in git
cd ../agentA_work
python3 enl.py            # enlargement table (the headline)
python3 exhaust.py 6      # exhaustive code supports of the region  -> supports.json
python3 exh6.py           # exact HNF over every <=6 violation set
python3 gmax2.py 1 120 <state.json>   # generic regional max-satisfy for any state
python3 regsolve2.py / zsolve.py / diag2.py / tgrow2.py   # region model, HNF, decompiler
```

## Established (my own computation, exact, no floating point)
1. **The residual region is an exactly-affine model.** `regsolve2.py` picks knobs so every
   region atom is affine in them (<=1 knob per monomial), so no equation is approximated.
2. At the 39,026 witness: region 33 atoms / 39 eqs, +a37887,+a41906 -> 35/41, 22 knobs.
   The whole affine system is **Q-CONSISTENT with rank = #knobs**, so it has a UNIQUE
   rational solution W, and W is non-integral in exactly 5 coordinates with denominators
   x642:2458959, x1329:p, x9413:p, x10903:p, x17325:p*2458959.
   => the region's full-solve conditions are exactly
   **7376877 | (x7068-x2099), p | x9118, p | x8731, p | x28730.**
3. **Enlarging the movable variables does NOT weaken the bound** (the question I was set):
   freeing x9118, x8731, x2099, x7068 (54 atoms / 20 knobs / 110 eqs) and the targeted
   growth to 55 knobs / 514 eqs both keep rank == #knobs and the SAME five denominators.
   In the enlarged region ISD (48k trials) finds **no code support of weight <= 6 at all**,
   so every integer knob vector violates >= 7 rows there.  Enlargement makes it stronger.
4. Exhaustive support enumeration of the 22-knob region (`exhaust.py`, information-set
   argument, provably complete): min support 5; 11,628 weight-5 and 27,303 weight-6
   supports, verified real over Q.  `exh6.py` tests all 193,971 admissible <=6 violation
   sets by mod-p filter + exact HNF.  [see runs/exh6.log for the verdict]
5. **New census**: 1,562 atoms carry a private variable, **326 with granularity 1** (fully
   free atom value).  The lab's s10/handles.py saw only 1,249 (free inputs only, all
   p-quantised).  None of the 326 lies in the twelve residual equations.
6. **mod9118_0 (39,009) decompiled exactly**: its entire residual is
   a21617 = c1*x14623 + c2*x27522 (mod p) and a29539 = 25692874*(x14853 - x1308) (mod p);
   x14623 and x14853 are FREE INPUTS and neither is in the other's cone (`cone2.py`).
   Shifting them zeroes both residues mod p but costs heavily downstream (38,975, exp1.py).

## Highest-value next experiment
Break C2 by moving x28730 off its residue mod p while holding a37887 = 0.  a37887 is a
CHECK in ONE equation (8680), it is quadratic in the knobs, and EIGHT region variables
(x950,x6947,x9629,x15120,x23754,x33168,x35619,x28730) have a37887 as their ONLY external
atom.  My linear model deliberately drops eq 8680; solving the affine 39-row system jointly
with the single quadratic a37887(delta) = 0 (Groebner/resultant in <=8 variables) is the one
place in this region where the lockstep rank argument does not apply.
