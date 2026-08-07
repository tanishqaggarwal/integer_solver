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


## Single highest-value next experiment
The 39,026 basin is closed.  Attack the CANONICAL basin with the same machinery: at
mod9118_0 the whole residual is a21617 = c1*x14623 + c2*x27522 (mod p) and
a29539 = 25692874*(x14853 - x1308) (mod p), with x14623 and x14853 FREE INPUTS and neither
in the other's cone.  Run `full31`-style exhaustive coset decoding on THAT region
(76 atoms / 65 knobs / 98 eqs, only two non-integral knobs x5040, x30163) instead of the
39,026 region — it is the only region found whose obstruction is 2 congruences rather than
2 congruences plus a rank-7 wall.

## ============ ADDENDUM ============

### Priority 1 — canonical basin at mod9118_0 (pure integer linear algebra)
Region: 89 nontrivial affine rows, 65 knobs, rank 65 = #knobs, Q-consistent with a UNIQUE
rational solution W that is non-integral in exactly two coordinates (x5040, x30163, both
denominator p).  Same lemma as the 39,026 theorem: every integer point's violated set D
must contain a code support.
* `mindist.py` — minimum support weight = minimum number of linearly DEPENDENT COLUMNS of
  the parity check H (24 x 89, = left kernel of N).  Exhaustive for |D| <= 3: **none**.
  Smallest single-knob column support is 10.
* `modpobs.py` — the decisive filter, mirroring what killed all 38,760 candidates in the
  39,026 region.  Reduce mod p: 33 of the 65 knob columns VANISH mod p (the p-quantised
  handles), rank(N mod p) = 32, left-kernel dim w = 57, and the syndrome g = Wb.B is
  NONZERO.  The retained rows V\D are mod-p consistent iff g lies in the span of the |D|
  columns of Wb indexed by D.  Exhaustive: **no |D| <= 3 works.**
  This reduces the whole question to one minimum-weight syndrome-decoding problem over
  F_p: min |supp(x)| with Wb*x = g, i.e. the coset minimum weight of a [89, 32] code.
* `canon2.py` — VERDICT: 10,635 information-set trials, **zero supports of weight <= 6**;
  lightest support seen has weight 10 (= the smallest single-knob column support).
  Detection probability 4.9% per trial for weight 6 => P(missed) <= e^-534.
* `prange.py` — the independent mod-p route: 400 Prange trials, 190 solvable, **none of
  weight <= 6**; detection probability 0.0624 per trial => P(missed) <= 4.8e-6.
  (Weights SEEN are upper bounds on the true minimum, not lower bounds; the sound
  inference is the absence of weight <= 6.)  Rigorous exhaustive floors: no dependent
  column pair/triple of H, and no mod-p-consistent |D| <= 4.
  **VERDICT: the canonical basin at mod9118_0 cannot beat 39,026 — closed.**

### Status of the three basins I have measured
| basin | rows/knobs | min violated (integer) | best score |
|---|---|---|---|
| 39,026 witness region | 31 / 22 | **7, proven exhaustively** | 39,026 |
| enlarged 39,026 region | 92 / 20 | >= 7 (no support of weight <= 6 in 48k trials) | 39,026 |
| canonical mod9118_0 | 89 / 65 | **10** (min support weight; P(miss) <= e^-534) | 39,009 |
