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


## ============ EQUATION-LEVEL ADDENDUM (the campaign's central question) ============
The shared hidden assumption was NOT "all atoms zero" in the rows — my rows were always
equation values — but in the REGION: the foreign atoms of every modelled equation were
required to stay zero.  The fix is the EQUATION closure A := atoms(eqs(A)), which admits
them as free cancellers.  `eqwin.py`, `eqwin2.py`, `eqbound.py`, `eqisd.py`, `boundary.py`,
`cancel.py`.

  lev  atoms  eqs  rows  knobs  rank  Qincons  violated   vars  excluded
   0     24    27    15     9      9     0        7        56     47
   2     88    94    68    32     32     0        7       202    170
   4    163   162   113    80     80     0        7       374    294
   6    235   230   163   109    109     0        7       537    428

* Every window is exactly affine — ZERO atoms nonlinear in the knobs at any level, so no
  equation is approximated and no knob is dropped for linearity.  All exclusions are
  variables touching an atom outside the window; that is the whole boundary.
* rank = #knobs and Q-consistency hold at every level => the uniqueness lemma applies and
  every integer point's violated set must contain an equation-level code support.
* L=2: EXHAUSTIVE mod-p floor >= 5; Prange 2025 solvable trials, none of weight <= 6,
  P(missed) <= 9.5e-208.  L=6: exhaustive floor >= 4; Prange 195 solvable trials, none of weight <= 6, P(missed) <= 9.8e-09.
* **>= 7 SURVIVES at equation level.**  Admitting 211 extra cancelling atoms and 100 extra
  knobs produces no cheaper point.
* `cancel.py`: 3,235 atoms occur in exactly one equation and **none carries a private
  handle** — no single-equation atom is independently settable anywhere in the instance,
  and none appears in any of the 7 failing equations (their atoms occur in 6-14 equations).
  The cheap-cancellation lever does not exist where it would matter.

### Jobs still tightening (results already established will only get sharper)
* `runs/eqisd6.log`  — independent condition (b), code-support ISD at L=6, 0.0981/trial.
* `runs/eqb2_exh5.log` — exhaustive mod-p floor |D| <= 5 at L=2 (raises >= 5 to >= 6).
* `runs/eqb6b.log`  — 4000 Prange trials at L=6 (drives P(missed) from 1e-8 towards 1e-180).
Re-run any of them with:  python3 eqbound.py <state> <level> 6 <trials> <exhaustive_k>


## ============ CORRECTION AND CEILING ============
* My claim "the excluded set shrinks as L grows" was **WRONG**.  It grows monotonically
  (47 -> 428 at L=6 -> 11,794 at L=344 -> 17,220 at L=791) and the excluded FRACTION stays
  near 0.80 throughout.  `fastgrow.py`, `runs/grow.log`, `growth.json`.
* **The exactly-affine regime ends at L=344.**  L=345 is the first level with an atom
  nonlinear in the knobs; by L=791 there are 1,287.  Closure (excluded = 0) needs the whole
  giant component (23,843 vars), of order 1,500+ levels.  The affine ceiling arrives first,
  so **the theorem cannot be made unconditional by raising L** — the conditionality is
  structural, not a compute budget.
* `eqisd.py` at L=6 was VACUOUS (0 of 3,300 uniform 109-subsets of 163 rows had rank 109).
  Use `eqmindist.py` instead: rigorous dependent-column floors + greedy information sets.
* Condition (a) final: L=2 exhaustive floor >= 5 and P(missed) <= 9.5e-208;
  L=6 P(missed) <= 6.7e-64; L=16 in `runs/eqb16.log`.

## Reconciliation with a global kernel computation (agent F's question, stated abstractly)
If a global computation finds a nontrivial kernel of the equation-atom incidence matrix M,
that is an ATOM vector b with M b = 0 — every equation satisfied by cancellation.  My
windows say no integer knob vector in them reaches fewer than 7 violated equations.  The
two are compatible and the reconciliation is forced: such a kernel vector is not in the
IMAGE of the knob-to-atom map.  My models make that concrete — the knob-to-atom map is
exactly affine with rank = #knobs, so its image is a single coset of a rank-|K| lattice in
Z^|A|, while ker(M) has dimension |A| - rank(M) which is far larger.  Realizability, not
cancellation, is what is missing.  Independent support: `cancel.py` shows 3,235 atoms occur
in exactly one equation and NONE carries a private handle, so the atoms a kernel vector
would need to move most freely are exactly the ones that cannot be moved independently.

## ============ SOUNDNESS OF THE EQUATION-LEVEL BOUND ============
Checked, not assumed (`soundness.py`, `soundness2.py`, `cancel.py`, `w6test.py`):
1. My code support ranges over **knob directions**, not atom vectors:
   C = {N u : u in Q^K}, K actual variables, N[e][j] = d(eq e)/d(knob j).  The knob->atom
   map has image a rank-<=|K| sublattice of Z^|A| — at L=6, 109 of 235 atoms are movable
   and 126 are FROZEN at 0; at L=16, 334 of 611 movable, 277 frozen.
2. The RAW relaxation over arbitrary atom vectors IS vacuous, confirmed independently in
   my parse: 3,235 atoms occur in exactly one equation, so min ||M a||_0 = 1.
3. A weight-1 support is not excluded by fiat, but does not occur: ZERO single-equation
   atoms lie in any of my windows, and globally NONE of the 3,235 single-equation atoms
   carries a private variable, so none is independently settable anywhere.  Rigorously the
   minimum support is >= 4 at L=6 (no 2- or 3-column dependency, 721,764 subsets) and 6 as
   observed.  The witness a39032 occurs in EIGHT equations in my parse and has no private
   variable.
4. Honest nuance: low-occupancy atoms DO drive my minimum support down — the observed 6 is
   exactly a35758's equation set (6 equations).  They do not drive it to 1.  All 582
   weight-6 supports fail the mod-p filter, which is what leaves 7.
Conclusion: the bound stands.  Realizability enters through the knob construction, which
is also the ker(M) reconciliation.

## Condition (b), final numbers
L=2 (68 rows, 32 knobs, H = 36x68): RIGOROUS no dependent column subset of size <= 4
  (814,385 exhausted) => minimum support >= 5; size-5 search running.
L=6 (163 rows, 109 knobs, H = 54x163): RIGOROUS no size <= 3 (721,764 exhausted) => >= 4;
  greedy information sets give 258/258 usable trials (vs 0/3,300 uniform), lightest
  support = 6, 582 weight-<=6 supports enumerated, ALL fail mod-p, none integral.
L=16 (486 rows, 334 knobs): EXHAUSTIVE mod-p no |D| <= 2 (117,855) => >= 3; Prange 250
  trials, lightest weight 7, P(missed) <= 1.6e-4; deeper campaigns in runs/eqb16b.log and
  runs/w6_16.log.
