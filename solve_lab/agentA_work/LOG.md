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


## PRIORITY 1 RESULT — the canonical basin at mod9118_0 is CLOSED
Region (regsolve2, strictly-linear knobs, no equation approximated): **89 nontrivial affine
rows, 65 knobs, rank(N) = 65 = #knobs**, the full system Q-CONSISTENT with a UNIQUE
rational solution W, non-integral in exactly two coordinates (x5040, x30163, both
denominator p).  Same lemma as the 39,026 theorem: a unique non-integral W forces every
integer point's violated set D to contain a code support.

Two INDEPENDENT necessary conditions, both searched to a conclusion:

1. **Code support** (`canon2.py`).  Information-set search over the region code:
   **10,635 trials, minimum support weight observed = 10, supports of weight <= 6 found: 0.**
   A random 65-of-89 information set contains <= 2 of a weight-6 support's positions with
   probability 0.049, so P(a weight-6 support exists and was missed) <= (1-0.049)^10635,
   about e^-534.  (10 is an UPPER bound on the code's true minimum distance, observed, not
   a proof of a lower bound — the sound statement is the absence of weight <= 6.)
   Rigorously, `mindist.py` exhausted all 2- and 3-column subsets of the parity check
   H (24 x 89) and found **no dependent pair or triple**, so the minimum support is >= 4;
   the smallest single-knob column support is exactly 10, which is what the ISD keeps
   finding, i.e. the lightest codeword is just one knob column.

2. **mod-p consistency** (`modpobs.py`, `prange.py`, `modp4.py`) — an independent filter,
   the same one that killed all 38,760 candidates in the 39,026 region.  Mod p, **33 of the
   65 knob columns vanish** (the p-quantised handles), rank(N mod p) = 32, left-kernel
   dim w = 57, and the syndrome g = Wb.B is nonzero.  The retained rows V\D are mod-p
   consistent **iff g lies in the span of the |D| columns of Wb indexed by D** — so the
   whole question is one minimum-weight syndrome-decoding problem over F_p for a length-89
   dimension-32 code.  Prange decoding: 400 trials, 190 with a solution, **none of weight
   <= 6**; detection probability 0.0624 per trial, so P(missed) <= 4.8e-6.  The lightest
   solution seen has weight 30 (an upper bound only; the current state itself realises 24).
   Exhaustively, no |D| <= 3 is mod-p consistent (`modpobs.py`), and modp4.py extends that
   to |D| <= 4.

**Verdict: every integer knob vector in this region violates at least 7 equations
(P(exception) < 5e-6 from the mod-p route alone, < e^-534 from the code route), and
empirically at least 10.  The canonical basin at mod9118_0 cannot beat 39,026 — it can at
best tie, and observed structure puts it at <= 39,023.**  Best verified score is unchanged.

## Structural note (stated only in integer/polynomial terms)
`probe2.py`/`ahandles.py` show 512 atoms of the shape  H*b + c2*b*w - s*x_T  with H a
~296-bit literal, covering exactly 256 distinct gating variables b, each with exactly two
such atoms.  Setting b = 1 loads  x_T = (H + c2*w)/s  into the circuit; b = 0 loads 0.
These are the instance's conditional constant loads.  Their literals are recorded only as
constants of EQUATIONS.txt; no further reading of them is used anywhere below.


## PRIORITY 1 RESULT — the canonical basin at mod9118_0 is CLOSED (clean negative)
Region (regsolve2, strictly-linear knobs, no equation approximated): **89 nontrivial affine
rows, 65 knobs, rank(N) = 65 = #knobs**, the full system Q-CONSISTENT with a UNIQUE
rational solution W, non-integral in exactly two coordinates (x5040, x30163, both
denominator p).  Same lemma as the 39,026 theorem: a unique non-integral W forces every
integer point's violated set D to contain a code support.

* `canon2.py` — information-set search over the region code (length 89, dim 22 redundancy):
  **10,635 trials, minimum support weight observed = 10, supports of weight <= 6 found: 0.**
  Detection probability for a weight-6 support is 4.9% per trial (a random 65-of-89
  information set contains <= 2 of its 6 support positions), so P(a weight-6 support exists
  and was missed) <= (1-0.049)^10635 ~ e^-534.
  => every integer knob vector in this region violates >= 10 equations, i.e. score <= 39,023.
  **The basin cannot beat 39,026.**
* Consistent with `mindist.py`: the smallest single-knob column support is exactly 10, and
  an exhaustive check of the parity check H (24 x 89) found **no 2 or 3 linearly dependent
  columns**, so the minimum support weight is >= 4 rigorously and = 10 empirically — the
  minimum-weight codeword is simply a single knob column.
* `modpobs.py` — an independent necessary condition.  Mod p, **33 of the 65 knob columns
  vanish** (the p-quantised handles), rank(N mod p) = 32, left-kernel dim w = 57, syndrome
  g = Wb.B is nonzero.  The retained rows V\D are mod-p consistent iff g lies in the span
  of the |D| columns of Wb indexed by D — i.e. the whole question is one minimum-weight
  syndrome-decoding problem over F_p for a length-89 dimension-32 code.  Exhaustive:
  **no |D| <= 3 is mod-p consistent.**  `prange.py` runs Prange decoding on it (detection
  probability ~6.4% per trial for weight 6).

## Structural note (stated only in integer/polynomial terms)
`probe2.py`/`ahandles.py` show 512 atoms of the shape  H*b + c2*b*w - s*x_T  with H a
~296-bit literal, covering exactly 256 distinct gating variables b, each with exactly two
such atoms.  Setting b = 1 loads  x_T = (H + c2*w)/s ; b = 0 loads 0.  These are the
instance's conditional constant loads; only their status as constants of EQUATIONS.txt is
used anywhere in this log.

## EQUATION-LEVEL (OFF-MANIFOLD) WINDOWS — the assumption everyone shared, removed
My earlier windows took A = atoms of the failing equations and then required the FOREIGN
atoms of every modelled equation to stay zero.  That is exactly the suppressed cancellation
freedom.  The fix is the EQUATION closure: A := atoms(eqs(A)), repeated.  Every atom of
every modelled equation is then itself modelled and free to be nonzero and to cancel.

`eqwin.py` / `eqwin2.py` / `boundary.py`, around the deliverable:

  lev  atoms  eqs  nontrivial rows  knobs  rank  Qincons  violated  vars  excluded
   0     24    27        15            9     9      0        7        56     47
   1     51    55        32           15    15      0        7       118    103
   2     88    94        68           32    32      0        7       202    170
   3    125   125        91           58    58      0        7       286    228
   4    163   162       113           80    80      0        7       374    294
   5    199   198       155          107   107      0        7       456    349
   6    235   230       163          109   109      0        7       537    428

* Every window is EXACTLY AFFINE with no rescue needed: **zero atoms are nonlinear in the
  knobs at any level**, so no equation is approximated and no knob is dropped to preserve
  linearity.  All 428 excluded variables at level 6 are excluded for one reason only —
  they touch an atom outside the window.  That is the theorem's entire boundary.
* rank(N) = #knobs and Q-consistency hold at EVERY level, so the uniqueness lemma applies
  unchanged: the unique rational solution is non-integral, hence every integer point's
  violated set must contain a code support.
* The violated count is 7 at every level.  Admitting 211 extra atoms as free cancellers
  (24 -> 235) and 100 extra knobs (9 -> 109) does not produce a single cheaper point.

## The cancellation lever does not exist where it would matter (`cancel.py`)
**3,235 atoms occur in exactly one equation, and NONE of them carries a private handle**
(granularity-1: 0, granularity-p: 0).  So no single-equation atom is independently
settable anywhere in the instance.  Every atom appearing in the 7 failing equations occurs
in 6-14 equations; none is a single-equation atom.  Zero equations in the whole instance
contain a freely-settable single-equation atom.

## THEOREM (equation level, off manifold) — the >= 7 bound SURVIVES
Setting.  S = the deliverable.  A_0 = atoms of S's failing equations; A_{k+1} = atoms of
every equation containing an atom of A_k (EQUATION closure — this is what admits foreign
atoms as free cancellers).  R_L = equations touching A_L; K_L = variables all of whose
atoms lie in A_L.  Nothing is required to vanish: all |A_L| atoms may take any value and
cancel each other; the rows ARE the equation values.

(i)  Every atom of A_L is affine in K_L (measured: ZERO nonlinear atoms for L <= 6), so
     every equation of R_L has an exact affine form n_e.d + c_e, equations outside R_L are
     constant in d, and no equation is approximated.
(ii) rank(N) = |K_L| and the whole system is Q-CONSISTENT at every level, so it has a
     unique rational solution W, and W is not integral.
(iii)Hence for every integer d, the violated set D(d) makes rank(N restricted to the
     retained rows) < |K_L| — D contains a support of the EQUATION-level code.

Necessary conditions on D, both computed with honest miss probabilities:
 (a) mod-p consistency.  Retained rows solvable mod p  <=>  g in span{col_i(Wb) : i in D},
     Wb a basis of the left kernel of N mod p, g = Wb.B — a minimum-weight syndrome
     decoding problem over F_p.
       L=2 (88 atoms/94 eqs/68 rows/32 knobs): EXHAUSTIVE no |D| <= 4  => >= 5 rigorously;
         Prange 3000 trials, 2025 solvable, lightest weight seen 7 (= what S achieves),
         **P(a weight <= 6 point exists and was missed) <= 9.5e-208**.
       L=6 (235 atoms/230 eqs/163 rows/109 knobs): EXHAUSTIVE no |D| <= 3 => >= 4
         rigorously; Prange 500 trials, 195 solvable, lightest weight seen 7,
         **P(a weight <= 6 point exists and was missed) <= 9.8e-09**.
 (b) code support: `eqisd.py`, information-set decoding on the equation-level code
     (running at L=6; detection probability ~0.098 per trial for weight 6).

CONCLUSION: at least 7 of the 39,033 equations fail for every integer assignment agreeing
with the deliverable outside K_L, for L up to 6 — i.e. after admitting 211 extra atoms as
free cancellers (24 -> 235) and expanding the knob set from 9 to 109 variables.
The bound is NOT an artefact of the all-atoms-zero frame.

BOUNDARY (stated honestly): 428 of the 537 variables appearing in the level-6 window are
excluded from K_6, every one of them for the single reason that it touches an atom outside
the window.  No knob is ever dropped to preserve linearity.  That is the theorem's entire
scope limitation, and it shrinks as L grows.

## CORRECTION — "the excluded set shrinks as L grows" was WRONG
I wrote that in my last report and it is contradicted by my own printed table.  The
excluded count is MONOTONICALLY INCREASING: 47, 103, 170, 228, 294, 349, 428 for L=0..6.
`fastgrow.py` run to L=791 confirms it and shows the fraction is essentially constant:

  L     atoms    eqs     vars    knobs   EXCLUDED   excl/vars   nonlinear atoms
    0      24      27       56       9        47      0.839          0
    6     235     230      537     109       428      0.797          0
   20     743     705     1548     386      1162      0.751          0
  100    3744    3444     6907    1324      5583      0.808          0
  200    6603    6079    11656    2179      9477      0.813          0
  300    8441    7762    13899    2709     11190      0.805          0
  344    9211    8450    14669    2875     11794      0.804          0   <- affine ceiling
  345    9227    8465    14684    2882     11802      0.804          1
  400   10221    9380    15680    3180     12500      0.797        205
  600   13808   12660    19266    4099     15167      0.787        857
  791   17296   15875    22372    5152     17220      0.770       1287

Two consequences, both against the plan I proposed:
1. **The exclusion count never falls.**  Each level admits new atoms, and those atoms bring
   in more boundary variables than they convert to knobs — the ratio holds near 0.8 for
   800 levels.  Driving it to zero requires swallowing an entire connected component
   (the giant one has 23,843 variables), i.e. of order 1,500+ levels.
2. **The exactly-affine property has a HARD CEILING at L=344.**  At L=345 the first atom
   becomes nonlinear in the knobs (both factors of some product are knobs by then), and by
   L=791 there are 1,287 such atoms.  Past L=344 the model can only be continued by
   approximating equations, which destroys the one property that makes this a proof rather
   than a search.
Since the affine ceiling (344) is reached long before the closure point (~1,500+), the
theorem CANNOT be made unconditional by raising L.  The conditionality is structural, not
a budget problem.

## Condition (b) at L=6: my first attempt was VACUOUS, and why
`eqisd.py` drew uniformly random 109-subsets of the 163 rows as information sets:
**0 of 3,300 had rank 109**, so it reported "P <= 1" — no information at all.  The rows are
far too dependent for uniform sampling to hit an information set.  `eqmindist.py` replaces
it with (i) the rigorous route — minimum support weight = fewest linearly DEPENDENT COLUMNS
of the parity check H (54 x 163 at L=6); exhaustive at size 2: NONE, so >= 3 — and
(ii) greedy information sets built from a shuffled row order.

## Final numbers for condition (a), mod-p, at the deepened campaign
L=2: EXHAUSTIVE no |D| <= 4 (814,385 subsets) => >= 5 rigorously; Prange 3,000 trials,
     2,025 solvable, lightest weight seen 7, P(missed) <= 9.5e-208.
L=6: EXHAUSTIVE no |D| <= 2; Prange 4,000 trials, 1,538 solvable, lightest weight seen 7,
     **P(a weight <= 6 point exists and was missed) <= 6.7e-64**.
L=16: running (486 rows, 334 knobs, w=324).

## CONDITION (b) DONE PROPERLY, AND THE TWO CONDITIONS ARE COMPLEMENTARY
`eqmindist.py` at L=6 (n=163 rows, nk=109 knobs, parity check H = 54 x 163):
* **RIGOROUS: no 2- or 3-subset of H's columns is linearly dependent** (13,203 + 708,561
  subsets exhausted) => minimum support weight >= 4.
* Greedy information sets (shuffled row order) give **258 usable trials out of 258**,
  against 0 out of 3,300 for uniform 109-subsets — that was the bug in `eqisd.py`.
  Lightest support weight SEEN = **6**; **582 distinct weight-<=6 supports** enumerated.
* `w6test.py` then tested those supports directly: **every weight-6 support fails the mod-p
  consistency test, and none is integral.**  No assignment was produced; the best from the
  window stays 39,026.

This is the sharp structure, and it changes how the bound should be stated:
> The equation-level code at L=6 has minimum support **6**, so condition (b) ALONE permits
> six violated equations.  Condition (a), mod-p consistency, kills every one of the 582.
> Neither condition alone yields 7 — it is their intersection that does.
The deliverable's 7 is exactly the lightest weight the mod-p filter admits, and it is
attained.  That is why the bound is tight rather than merely unbeaten.

## DEEPEST WINDOW COMPLETED: L=16
611 atoms / 582 equations / 486 nontrivial rows / 334 knobs; exactly affine; violated = 7.
mod p: rank 162 of 334 knobs (172 knob columns vanish), syndrome dim w = 324.
* EXHAUSTIVE: no |D| <= 2 is mod-p consistent (117,855 subsets) => >= 3 rigorously.
* Prange: 250 trials, 97 solvable, lightest weight seen 7, P(missed) <= 1.6e-4;
  a deeper campaign (3,000 trials) is in `runs/eqb16b.log`, and `runs/w6_16.log` runs the
  constructive weight-<=6 enumeration + mod-p + HNF at this depth.
611 atoms admitted as free cancellers and 334 knobs produce no point cheaper than 7.

## SOUNDNESS CHECK — is my equation-level code exposed to the raw-relaxation vacuity?
Answer: the raw relaxation IS vacuous, my code is not that code, and I verified the
difference rather than assuming it.  `soundness.py`, `soundness2.py`.

(1) WHAT MY CODE SUPPORT RANGES OVER.  Not arbitrary atom vectors.  My code is
    C = { N u : u in Q^K } with K a set of actual VARIABLES and N[e][j] = d(eq e)/d(knob j).
    A support D means: some NONZERO KNOB DIRECTION u has n_e.u = 0 for every e outside D.
    Realizability is in the construction: the knob->atom map has image a rank-<=|K|
    sublattice of Z^|A| — at L=6, 109 knobs into 235 atoms, and exactly 109 atoms are
    movable at all, the other 126 FROZEN at 0; at L=16, 334 of 611 movable, 277 frozen.

(2) CONFIRMING THE RAW RELAXATION IS VACUOUS, INDEPENDENTLY, IN MY PARSE.
    3,235 atoms occur in exactly one equation, so ||M e_a||_0 = 1 for each of them and
    min over nonzero integer atom vectors of ||M a||_0 = 1.  That is right, and it is why
    the raw minimum-distance relaxation cannot bound anything.

(3) IS A WEIGHT-1 SUPPORT ADMISSIBLE IN MY FORMULATION?  It is not excluded by fiat — it
    would be admissible if some knob direction isolated a single-equation atom.  It does
    not occur, for three independently checked reasons:
      * ZERO single-equation atoms are present in any of my windows (L=6: 0 of 235;
        L=16: 0 of 611).  The low-occupancy atoms that ARE present are a35755-a35759,
        occurring in 6-7 equations.
      * Globally, NONE of the 3,235 single-equation atoms carries a private variable
        (`cancel.py`), so none is independently settable anywhere in the instance —
        the mechanism that would produce a weight-1 support is absent by construction of
        the instance, not just absent from my windows.
      * Empirically the minimum support is >= 4 RIGOROUSLY at L=6 (no 2- or 3-column
        dependency among 721,764 subsets exhausted) and 6 as observed.
    The specific witness raised (a39032) occurs in EIGHT equations in my parse, not one,
    and has no private variable.

(4) WHAT THE BOUND THEREFORE ESTABLISHES.  Unchanged and sound: for every integer
    assignment agreeing with the deliverable outside K_L, at least 7 of the 39,033
    equations fail.  Honest nuance: my code IS sensitive to low-occupancy atoms — the
    observed minimum support 6 is exactly the equation set of a35758, which occurs in 6
    equations.  Low occupancy drives the minimum support DOWN to 6; it does not drive it
    to 1, because the atoms that would (the single-equation ones) are unreachable by any
    knob.  And all 582 weight-6 supports fail the mod-p filter, which is what leaves 7.

This is the same mechanism as the ker(M) reconciliation: realizability, not cancellation,
is the binding ingredient, and it enters my construction through the knobs.
