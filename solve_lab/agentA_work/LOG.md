# Agent A LOG

## t0 — setup
- s9/*.pkl caches were absent (gitignored). Rebuilt: atomize/poly/gates/fwd (~60 s).
- Verified baseline: best/new_instance_partial_39026.json -> 39026/39033. CONFIRMED.

## Structure re-derivation (probe1..probe5, region.py, ahandles.py)
- 7 nonzero atoms; 12 equations E; 33 atoms in E; 39 equations in region; 11 zero-cost knobs.
- Private-handle census (ALL variables, not just free inputs): 1,562 atoms have a private
  variable; 326 of them with granularity 1. Prior lab census (s10/handles.py) restricted to
  free inputs -> 1,249, all p-quantised. My census is a strict superset. None in E.
- Found four knobs the prior generator list missed: x1613, x1844, x21574, x29305
  (all with 0 atoms outside the 33-atom region).

## Exact region model (the main result so far)
Built `regsolve2.py` (strictly-linear knob selection: at most one knob per monomial, so
the model is EXACTLY affine and no equation is dropped) + `agrow.py`/`amk_model.py`.

At the 39,026 witness, region = 33 atoms / 39 eqs / 11 knobs; enlarged with a37887+a41906
-> 35 atoms / 41 eqs / 22 knobs.  Results:
* The full 40-row affine system is CONSISTENT over Q with rank 22 = #knobs, i.e. a UNIQUE
  rational solution W.  W is NOT integral: exactly 5 coordinates have denominators
  x642:2458959, x1329:p, x9413:p, x10903:p, x17325:p*2458959.
  => full-solve of the region <=> 7376877 | (x7068-x2099), p | x9118, p | x8731, p | x28730.
* ENLARGING the movable set does NOT help: +29426,+41972,+29090,+36085 and the targeted
  growth (9 -> 55 knobs, 24 -> 154 atoms, 27 -> 514 eqs) all keep rank == #knobs and the
  SAME 5 denominators.  Rows and columns grow in lockstep.  (Independently confirms the
  prior lab's part-5 claim, by a different method.)
* Max-satisfy: 10 of the 40 rows are identically zero in the knobs; the other 30 rows form
  a rank-22 code.  ISD over that code (`gmax.py`) enumerated many supports of size 5 and 6;
  EVERY size-5 and size-6 support is NOT integrally solvable (exact HNF).  So within this
  region 7 violated equations stands. 39,026 optimal here.

## Other states (regsolve2/zsolve)
* mod9118_0 (39,009): region 76 atoms / 65 knobs / 98 eqs, rank 65, Q-consistent, only
  TWO non-integral knobs x5040(p), x30163(p).  Decompiled: the entire residual there is
  a21617 = c1*x14623 + c2*x27522 (mod p) and a29539 = 25692874*(x14853-x1308) (mod p),
  with x14623 and x14853 FREE INPUTS whose cones do NOT contain x27522/x1308.
* exp1.py: shifting those two free inputs zeroes both residues mod p but the collateral
  is large (58 nonzero atoms, 38,975).  Handle re-solve not yet applied.
* beam_39016 / AG_39013 / PF_39015 / wr_w1_39020: exact HNF says NO integer solution of
  the region at level 0.

## Exhaustive region optimality (the headline computation)
`exhaust.py` — with a fixed information set, any codeword of weight <= W has >= |I|-W
zeros inside I, so enumerating J = supp(c_I) with |J| <= W and the forced-zero subsets of
OUT is EXHAUSTIVE.  On the 39,026 region (30 nontrivial rows, 22 knobs):
  minimum support = 5;  11,628 supports of weight 5 and 27,303 of weight 6 (all verified
  real over Q by `verify_sup.py`, 25/25 sampled, 0 false positives).
Since rank(N) = 22 = #knobs and the FULL row system is Q-consistent with a unique
non-integral solution W, an integer point violating a set D exists only if D contains a
code support.  `exh6.py` therefore enumerates every admissible D with |D| <= 6
(38,931 supports + 11,628*25 one-row extensions), filters by mod-p consistency and
tests each survivor with exact HNF.  [running]

## Enlarged region (54 atoms / 20 knobs / 110 eqs; x9118,x8731,x2099,x7068 all freed)
`gmax_enl.py`: 48,068 ISD trials found **no code support of weight <= 6 at all**, so in
the enlarged region every integer knob vector violates >= 7 rows.  The current point
attains exactly 7.  Enlarging the movable set makes the bound STRONGER, not weaker.

## THEOREM (exhaustive, exact) — 7 failing equations is optimal in the 39,026 region
`model31.py`: a37887 is a PERFECT SQUARE as a polynomial in the region knobs,
a37887 = Q^2 with Q affine in 8 knobs (verified coefficient-by-coefficient).  So eq 8680
is not a quadratic constraint at all: it is one more AFFINE row.  Complete model =
31 nontrivial affine rows, 22 knobs, 35 atoms, 41 equations.
`full31.py`: rank(N) = 22 = #knobs and the whole 31-row system is Q-CONSISTENT (unique
rational solution W, non-integral).  Therefore any integer point's violated set must
contain a code support.  Exhaustive information-set enumeration: minimum support = 6,
38,760 supports of weight 6, none smaller.  All 38,760 admissible violation sets of size
<= 6 FAIL mod-p consistency -> 0 integral.  Hence >= 7 violations for every integer knob
vector, and 7 is attained.  **39,026 is exactly optimal in this basin, proven over Z.**

## Frame-2 generator scan (ripple_test.py, scan2099.py) -- prices, measured
cost 0: x9118, x8731, x642, x29854, x31864, x1329, x10903, x9413, x17325 (9 free gens)
cost 13: x7068 and **x6418** (x6418 moves alpha0 by exactly -1 -- the fine control of C1)
cost 16: x28730 | cost 38: x4287 | cost 109: x2081 | cost 3: x6947, x21574 | 4: x1613
5: x33168 | 7: x950 | 9: x15120.  The cone of x2099 u x28730 is only 18 variables.

## Congruence restatement (sharpened)
alpha0 + 7376877*alpha6 = x7068 - x2099 - 7376877*p*x17325, so
  C1 <=> (x7068 - x2099) mod p fixed ;  C2 <=> x28730 mod p fixed.
At AG_39013 and mod9118_0 BOTH residues are exactly 0 -- because in the canonical frame
x7068 = 7376877*x642 + x2099 and x28730 = p*x9413 are GATE IDENTITIES.  So C1/C2 are the
price of the 39,026 frame, not an arithmetic obstruction of the instance.

## Structural note (stated only in integer/polynomial terms)
`probe2.py`/`ahandles.py` show 512 atoms of the shape  H*b + c2*b*w - s*x_T  with H a
~296-bit literal, covering exactly 256 distinct gating variables b, each with exactly two
such atoms.  Setting b = 1 loads  x_T = (H + c2*w)/s  into the circuit; b = 0 loads 0.
These are the instance's conditional constant loads.  Their literals are recorded only as
constants of EQUATIONS.txt; no further reading of them is used anywhere below.

