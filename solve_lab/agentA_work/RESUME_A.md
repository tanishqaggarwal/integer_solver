# Agent A — RESUME (exact integer linear algebra / lattice)   [FINAL]

## Best verified: 39,026 — solve_lab/agentA_work/A_best_39026.json
Identical to the lab baseline; re-verified with solve_lab/checker.py (39026/39033,
failing [12231,12270,12350,14584,18673,22044,29125]).  I did not beat it.

## Re-entry (the s9/*.pkl caches are gitignored and were MISSING — rebuild first)
cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py
cd ../agentA_work && python3 enl.py ; python3 model31.py ; python3 full31.py 6

## THEOREM (exhaustive, exact, no floating point)
Region of the 39,026 residual = 35 atoms / 41 equations / 22 knobs, chosen so every atom is
AFFINE in the knobs (regsolve2.pick_knobs: <=1 knob per monomial).  a37887 is a perfect
square, a37887 = Q^2 with Q affine, so eq 8680 is one more affine row -> 31 nontrivial rows.
  rank(N) = 22 = #knobs and the FULL 31-row system is Q-CONSISTENT: a UNIQUE rational
  solution W, non-integral in exactly 5 coordinates (denominators 2458959, p, p, p,
  p*2458959).  Hence every integer point's violated set must contain a code support.
  Exhaustive information-set enumeration: minimum support = 6, 38,760 supports of weight 6,
  none smaller.  All 38,760 admissible <=6 violation sets fail mod-p consistency.
  ==> >= 7 violated for every integer knob vector; 7 attained.  39,026 EXACTLY OPTIMAL here.
Enlarging the movable set does not weaken it: freeing x9118, x8731, x2099, x7068
(54 atoms / 20 knobs / 110 eqs) and targeted growth to 55 knobs / 514 eqs keep
rank == #knobs and the SAME five denominators; in the enlarged region ISD (48k trials)
finds no support of weight <= 6 at all.

## The obstruction, exactly
Reachable alpha (the seven atom values) = { alpha1 = K2 (mod p) } and
{ alpha0 + 7376877*alpha6 = C0 (mod p) }; 12 rows of rank 7 -> all 12 need alpha = 0 mod p.
Breaking one congruence buys exactly 1 equation (39,027); breaking both gives 39,033.
Sharpened:  C1 <=> (x7068 - x2099) mod p ;  C2 <=> x28730 mod p.
Both are GATE IDENTITIES of the canonical frame (x7068 = 7376877*x642 + x2099,
x28730 = p*x9413) and are exactly 0 at AG_39013/mod9118_0 — so they are the price of the
39,026 frame, not an arithmetic obstruction of the instance.
Full-solve conditions of the region: 7376877 | (x7068-x2099), p | x9118, p | x8731, p | x28730.
Measured levers (frame-2 ripple): x6418 = 13 (moves alpha0 by exactly -1), x7068 = 13,
x28730 = 16, x4287 = 38, x2081 = 109; nine zero-cost generators
(x9118,x8731,x642,x29854,x31864,x1329,x10903,x9413,x17325).

## Curve findings (independent)
b1 = 16469404786402603598127746642812631771238817117136746083575784224822817945026 from the
pinned point (x12186,x16742); b1/7 square but not cube -> cubic twist, NOT secp256k1;
group order N' = 109903*12977017*383229727*211853322379233867315890044223858703031485253961775684523,
COMPOSITE.  Over all 2817^2 literal pairs the max multiplicity of b = Y^2-X^3 is 1;
zero literal pairs on y^2=x^3+7 or on x-shifted versions.  Measured addition offset
K = 9941218437270274411588837402253980960504855302801171729868401674372857777188.

## Single highest-value next experiment
The 39,026 basin is closed.  Attack the CANONICAL basin with the same machinery: at
mod9118_0 the whole residual is a21617 = c1*x14623 + c2*x27522 (mod p) and
a29539 = 25692874*(x14853 - x1308) (mod p), with x14623 and x14853 FREE INPUTS and neither
in the other's cone.  Run `full31`-style exhaustive coset decoding on THAT region
(76 atoms / 65 knobs / 98 eqs, only two non-integral knobs x5040, x30163) instead of the
39,026 region — it is the only region found whose obstruction is 2 congruences rather than
2 congruences plus a rank-7 wall.

## ============ ADDENDUM (post-retasking) ============

### RETRACTION: my curve dissent was WRONG.  The instance IS a secp256k1 ECDLP.
My error: I tested the DEPRESSED short form y^2 = x^3 + b and never fitted the GENERAL
Weierstrass form (the a2*x^2 + a4*x terms); and I used arbitrary literal PAIRS instead of
the structurally-gated pin pairs.  Corrected derivation, entirely from my own model:
* `pins.py` — 512 load-pin atoms  H*bit + c2*bit*partner + s*x_target  covering exactly
  **256 gating bits, each with exactly 2 pins**.
* `leaf2.py` — the leaf coordinate is **-H mod p** (raw literal negated, NOT divided by s).
  222/256 inliers on one general Weierstrass curve (chance 3).  `weier.py` refit over all
  222: rank 3, **0 inconsistent rows**.
      a2 = 97553848499418123410591666447050222001188385549510401465815187079080512838891
      a4 = 114170008767671698752186727197936107864370654164657728518655355473804451402762
      a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490
* Depression X = x + a2/3:  **A = 0 exactly**,
      B = 64019533680030876408443198762210829058751700634554282185987325820393598524794,
  **B/7 is a sixth power mod p** => F_p-isomorphic to secp256k1; **n_secp*Pt = O for all
  222 points**, n_secp prime.  My earlier composite order N' was the curve you get by
  wrongly forcing a2 = a4 = 0 — a different curve.
* `chain2.py` — **all 256 leaf pins are exactly 2^i * G, i = 0..255, complete, no gaps**
  (isomorphism u = 12830242018875522506555146473674089970775060590290859819641972374662130570109;
  secp256k1's generator G is itself a leaf point).  189/222 points have their double in the set.
=> 256-bit double-and-add scalar multiplication on secp256k1; the selector bits ARE the
   scalar.  Agents C/D/G were right.  Artifacts: pins.json, weier.json, curve_final.json.

### Priority 1 (canonical basin at mod9118_0) — running
`canon2.py`: region 89 nontrivial rows / 65 knobs, rank 65 = #knobs, Q-consistent with a
unique non-integral W (only x5040, x30163, both denominator p).  So every integer point's
violated set must contain a code support; ISD searches for one of weight <= 6 and HNF-tests
each.  Log: runs/canon2.log.  (My earlier canon.py brute-forced 190,050 targeted drop-sets
without a filter and was too slow; canon2 replaces it.)
