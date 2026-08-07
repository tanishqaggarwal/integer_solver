# Agent H log — integer/polynomial analysis of EQUATIONS.txt

## Step 0
Read PROMPT.txt, RESUME.md, STATE.json.  Verified the claimed best:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> satisfied 39026/39033, failing [12231,12270,12350,14584,18673,22044,29125].  CONFIRMED.

## Step 1 — the frame (kept; this is the durable asset)
Independent parse (`model.py`): 42,267 atoms / 39,033 equations.  Orienting every atom of
syntactic form `x_t - rest` as the definition `x_t := rest` gives an ACYCLIC gate DAG over all
38,748 variables: 8,747 free inputs, 30,001 defined, 12,266 check atoms, ZERO cycles.
(The lab's prior frame had 7,273 free inputs and 1,800 variables inside gate cycles.)
Forward evaluation from ALL-ZERO free inputs scores **39,005** with only **5** nonzero check
atoms — checker-verified on `allzero_fwd.json`.

## Step 2 — decomposition hypothesis, REFUTED by my own measurement
eq-var bipartite graph: 1 component.  atom-eq graph: 1 giant + 3,234 singletons.
Free-input hypergraph (hyperedge = equation): 1 component, 8,747 vars / 18,248 equations;
20,785 equations are identically satisfied by forward evaluation.
Residual closure = 6,007 free inputs / 6,026 checks / 9,244 equations.
=> NO separator, NO block decomposition.  My assigned angle's central hypothesis is FALSE.

## Step 3 — the residual, in raw terms
Exactly three conditions remain:
  (i)  x_9274 = 1   (x_9274 = (x_7715 + x_34554) - x_7715*x_34554, over 256 selector variables)
  (ii) x_37892 = C1 (mod p)
  (iii)x_13682 = C2 (mod p)
with p = 115792089237316195423570985008687907853269984665640564039457584007908834671663,
C1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626,
C2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002.
Forced because a688 = 8863713*(x_18956 - C1) - x_14257 with x_14257 = p*x_7497 and gcd(8863713,p)=1;
a1618 = x_24468 - C2 - p*x_11436; a30982 = x_18956 - x_37892 - p*x_22820;
a30980 = x_24468 - x_13682 - 12354891*p*x_14393.

## Step 4 — how the 512 load-pin constants enter (measured, not assumed)
Each selector variable b carries exactly two pin atoms `b*(x_F - H) - c*x_handle` (512 constants,
`leafpins.json`).  They are combined by a tree of mux atoms.  For two active selectors the
combined pair (X3,Y3) is forced by the polynomial identities
   (X1+X2+X3+K)*(X2-X1)^2 = (Y2-Y1)^2 ,  (Y3+Y1)*(X2-X1) = (Y2-Y1)*(X1-X3)   (mod p)
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891.
Measured 56/56 on (u,w) selector pairs, with a 1-off control failing 56/56.
Over Z the triple additionally needs 2264251 | (x_15286/p), 2264251 = 11*43*4787; solving that
congruence for one lift parameter closes all three checks EXACTLY over Z (3/3).
**Consequence: the constants combine NONLINEARLY.  No linear relation among them expresses the
residual, so linear/lattice relation-finding over the constant table cannot be the operative tool.**

## Step 5 — tools built (kept)
`fast.py` incremental exact evaluator, 1.4 ms per move.
`close2.py` bottom-up cascade closer: closes the entire pin tree in ~32 exact assignments with NO
search (smallest-atom-support first, freeze each assigned variable).
`close3.py` adds an indirect phase.  `scan1.py`: every single-selector state closes to **39,018**
with exactly 4 nonzero atoms {30980, 30982, 36185, 40812}.
Beam search over repairs (`beam2.py`) plateaus at 39,012 and oscillates — the constructive closer
strictly dominates it.

## Step 6 — reframed task: integer relations over the constant table (all NEGATIVE, measured)
512 load-pin constants (`leafpins.json`), all distinct mod p, gcd 1.
- none equals C1 or C2 mod p; none equals C1±K or C2±K mod p
- all 130,816 pairwise differences: none equals C1 or C2 mod p
- all 130,816 pairwise sums:       none equals C1 or C2 mod p
- the additive-subset-sum formulation would have density 256/log2(p) = 1.000, the regime where
  lattice/low-density relation-finding provably does not apply — AND it is the wrong formulation
  anyway: the measured combination law is the NONLINEAR pair of polynomial identities in Step 4.

## Step 7 — frame B: the 39,026 witness reconstructed, then attacked with my extra knobs
`import26.py`: forward-evaluating the witness's free-input values in my frame gives 39,020 and
differs from the witness in EXACTLY 4 variables: x_642, x_28730, x_29854, x_31864.
`frameB.py`: detaching those 4 (8,751 free inputs, 12,270 checks) reproduces the witness EXACTLY —
score 39,026, 7 nonzero atoms {22229,22230,35758..35762}, 7 failing, 0 variables differing.
`climb.py`: exhaustive single-move scan over ALL 8,751 free inputs x {±1,±2,±p,±2p} = 70,008 moves,
plus every exact atom-solving root over the 7 nonzero atoms' supports.  **NO improvement.**
This covers the 1,478 free inputs the lab's prior frame did not have.
`comp.py`: region = 12 equations, 5 satisfied, 7 failing; exactly ONE atom (22231) has its whole
equation footprint inside the region and NO zero-collateral free input moves it; exactly NINE
zero-collateral knobs exist — {642, 1329, 8731, 9118, 9413, 10903, 17325, 29854, 31864}.
Exhaustive pair scan over those knobs x {±1,±p}: 576 pairs, best still 7 failing.
**This independently reproduces the lab's region census from a completely different orientation.**

## Step 8 — placement census in my own (non-witness) frame
`place.py`: all four ways of assigning the two congruence defects to atoms give
region == failing exactly (no cancellation at all): 15, 16, 22, 22 -> best 39,018.
`region.py`: that region admits ZERO compensator atoms.  My own frame's placement is strictly
worse than the witness's, and the witness's is the one my scans confirm at 7.

## Step 9 — THE BACKWARDS CLOSER (target region chosen first, values constructed to order)
`backreg.py` — global census of every target region generated by a single atom's equation
footprint: 23,059 distinct regions over the 23,844 atoms living in <=12 equations.
The DELIVERABLE's region measures |R|=12, |S|=8 (the seven plus a22231), **0 un-cancellable rows**.
Regions with better balance exist elsewhere in the instance (e.g. |R|=9,|S|=7,unc=0), but they
contain no defect-carrying atom, so they cannot hold the two congruences.

`backgrow.py` — greedy region growth seeded on the defect-carrying atoms, minimising |R|-|S|:
  witness seed  |R|=12 |S|= 8 balance 4 unc 0   ->  grows to |R|=23 |S|=21 balance 2 unc 0
  agent-H seed  |R|=25 |S|=17 balance 8 unc 3   ->  grows to |R|=27 |S|=21 balance 6
Combinatorial floor from the witness route: balance 2 + c 2 = **4 failing (score 39,029)**.

`backreal.py` — REALIZABILITY of that balance-2 region, measured in frame B.
29 free inputs reach S; exactly **9** move it with zero collateral, and they are the same nine
{642,1329,8731,9118,9413,10903,17325,29854,31864}.  Their images span only the ORIGINAL seven
atoms {22229,22230,35758..35762}.  The 14 extra atoms of the balance-2 region are unreachable.

`backopt.py` — exact integer optimum over that realizable lattice.  The 23 region rows are
JOINTLY LINEAR in the 9 knob offsets (verified).  Fraction-free integer elimination over every
row subset of size <= 9: **maximum 16 of 23 rows zeroable, i.e. 7 failing.  The witness already
attains it.**  Score 39,026.  Nothing was written above it.

`backprice.py` — price of each of the 14 unreachable atoms (best single move over its whole
support, deltas +-1, +-p):
  a22231 -> 12 failing | a22233 -> 14 | a22235 -> 15 | a19088, a19092 -> 17 | a10936, a19090 -> 19
  a10935, a10937, a19087, a19089, a19091, a22232, a22234 -> no move exists at all (empty support)
**None reduces the failing count below 7.**  (These are single-move upper bounds on cost; the
lab's deeper re-closure prices a22231 at 1 equation outside, and it buys exactly 1 — same verdict.)

### Verdict of the backwards experiment
Construction reaches the same 7 as search, and now explains why: the balance-2 region is a real
combinatorial object with 0 un-cancellable rows, but its atom-value lattice is 7-dimensional, not
21-dimensional.  The binding constraint is not the region's shape and not the atom count — it is
that only nine free inputs move the region with zero collateral, and their span is exactly the
seven atoms the deliverable already uses.  **|R|-|S| is the wrong objective; rank of the
realizable knob image is the right one, and it is 7 for every region reachable from the defect.**

## Step 10 — RANK-RAISING SWEEP (attack the knob image, not the region)
**Crossover argument (this is what makes the sweep exhaustive).** Adding one knob direction raises
the lattice rank by at most 1, so the number of zeroable rows rises by at most 1.  Every NEW
equation dragged into the region adds at least 1 to |R'|.  Net failing change >= (new eqs) - 1.
So a knob can pay ONLY if it drags in ZERO new equations.

`knobcensus.py` — collateral census over ALL 8,751 free inputs of frame B, deltas {1, p}:
    +0 new eqs :   9 knobs  {642,1329,8731,9118,9413,10903,17325,29854,31864}  (the base lattice)
    +1 new eq  :   1 knob   x_28730
    +2 new eqs :   1 knob   x_21574
    +3 / +4 / +5 : 1 / 1 / 2 knobs
    +7 and up  : thousands
**The +0 class is exactly the 9 knobs already in the lattice.  The sweep terminates at budget 1.**

Only 5 free inputs move a22231 (the sole compensator inside the region):
    x_28730 (+1 new eq), x_12553 (+15), x_4432 (+17), x_4287 (+36), x_2081 (+134).

`kernel22231.py` — could a COMBINATION cancel the collateral?  453 candidate knobs, 888 outside
atoms; over the 148 rows I could linearise exactly, rank goes 135 -> 136, i.e. a22231 looked
independent.  740 rows had no exact linear root, so that test was inconclusive by itself.
`combo5.py` settles it by direct evaluation (ground truth, no linearity assumption): all 16,806
integer combinations of the five movers with coefficients in [-3,3], plus a p-scaled pair sweep.
**Zero-collateral combinations found: 0.  Best failing: 7.**

`q37887.py` — the final reduction.  a22231 lives in 10 equations, ALL inside the region.
a37887 lives in **exactly one equation, eq 8680**.  Scanning all 8,751 free inputs: 17 move a37887;
exactly ONE of them (x_28730) moves nothing else outside the region — and it also moves a22231.
**There is no free input anywhere in the instance that moves a37887 without moving a22231**, so the
two cannot be made to cancel.  a22231 buys at most one row and costs exactly eq 8680.

### Verdict
The rank-7 realizable knob image cannot be raised at a profit.  The unique candidate is x_28730,
which raises the rank to 8 and drags in exactly one equation, eq 8680, and the trade is exactly
1-for-1.  Reached independently and from the free-input side, this lands on the same single row
another agent localised from the equation side.
**Scope of the claim:** this is exact for the witness placement, for single-direction augmentation
of the lattice, and for collateral measured at deltas {+-1, +-p}.  It is NOT a proof that no frame
anywhere admits 39,027; it is a proof that this region's knob image admits nothing better.

## Step 11 — PER-PLACEMENT SWEEP OF THE WHOLE CASCADE
`chain.py` extracts the cascade's pin chain from a one-selector closure: **20 pin atoms**, each with
its repairing free input.  `sweep.py` then puts the defect at each pin in turn (close everything
EXCEPT that atom) and prices the placement with the five-stage pipeline:
region R, inside-atoms S, balance deficit |R|-|S|, zero-collateral knob count, RANK of the
realizable knob image on the region rows, exact integer optimum (rows zeroable), failing, score.

carrier    nz  |R|  |S|  deficit  knobs  rank  zeroable  failing  score
WITNESS     7   12    8      4       9     7       5        7    39026   <-- rank > deficit
a688        4   15    6      9       2     2       0       15    39018
a23824      4   15    6      9       2     2       0       15    39018
a23826      4   15    6      9       2     2       0       15    39018
a14061      4   15    6      9       2     2       0       15    39018
a14063      4   15    6      9       2     2       0       15    39018
a12599      4   15    6      9       2     2       0       15    39018
a15100      4   15    6      9       2     2       0       15    39018
a9193       4   15    6      9       2     2       0       15    39018
a9195       4   15    6      9       2     2       0       15    39018
a6285       4   15    6      9       2     2       0       15    39018
a31928      4   15    6      9       2     2       0       15    39018
a26729      5   28   12     16       8     5      13       15    39018
a26731      5   28   14     14       9     6      13       15    39018
a1618       4   22   11     11       5     5       0       22    39011
a35267      5   26    7     19       3     3       0       26    39007
a31930      5   27    7     20       3     3       0       27    39006
a6283       5   28   11     17       3     3       0       28    39005
a21853      5   30    8     22       3     3       0       30    39003
a37733      5   30    7     23       3     3       0       30    39003
a37735      5   30    7     23       3     3       0       30    39003

**RESULT: the witness placement is the ONLY one in the entire cascade with knob-image rank
exceeding its balance deficit (7 > 4).  Every one of the 20 cascade pins has rank <= 6 against a
deficit of 9-23, and none scores above 39,018.  By my own criterion, 7 is the floor for the whole
cascade, not merely for one region.**

Structural reading: the cascade pins all sit in regions of 15-30 equations with only 2-3
zero-collateral knobs, because a pin's free input is consumed closing the pin itself.  The witness
placement is special because it sits where NINE free inputs act with zero collateral - the tail of
the tree, where the handles accumulate - and that is the only place in the instance where the
realizable image is wide enough to out-run the region it lives in.

**Scope:** exact for defect placement at a single cascade pin, single-selector closure, deltas
{+-1, +-p}, subset enumeration capped at min(knobs, rank, 6).  Not a proof that no frame anywhere
admits 39,027.

## Step 12 — HANDLE-CARRIER SWEEP (complete, not sampled)
**Correction to my own step-11 next-experiment note:** x_7497, x_11436, x_22820, x_14393 are NOT
solo handles — they occur in 2-3 check atoms each.  The real class is the free inputs occurring in
exactly ONE check atom: 1,865 of them, of which `handles.py` measures **1,143 with granularity
exactly p** and 722 dormant (zero effect).  Zero have any other granularity.

Construction (a handle carries the defect in multiples of p, not as a residue): from the closed
one-selector state, move the handle one step (its carrier atom shifts by m*p), FREEZE the handle so
it cannot absorb, re-close everything else, then price the placement.

`hsweep.py` — **all 1,147 carriers priced (1,143 solo + the 4 named non-solo), 122 s.  No sampling.**

  signature (nz,|R|,|S|,deficit,knobs,rank,failing,score)        count
  (3, 16,  4, 12, 2, 2, 16, 39017)                                440
  (4, 29,  5, 24, 3, 3, 29, 39004)                                 80
  (4, 30,  5, 25, 3, 3, 30, 39003)                                 56
  (4, 28,  5, 23, 3, 3, 28, 39005)                                 47
  (4, 27,  5, 22, 3, 3, 27, 39006)                                 37
  ... 39,017 is the maximum; the tail runs down to 38,997

  rank    : {2: 441, 3: 619, 4: 79, 5: 8}      max rank 5
  deficit : min 12, median ~24, max 31
  **carriers with rank > deficit: 0.  Minimum gap deficit - rank = 10.**

Stage B (exact integer optimum) was not needed: the winning criterion is rank > deficit and no
handle carrier comes within 10 of it, so none can beat 7 regardless of the optimum's exact value.
Nothing was written to disk; verifyE.py was not needed (no candidate state).

### Why handles are structurally worse than pins
A solo handle appears in exactly one check atom by definition, so breaking it yields a region with
only 4-6 inside-atoms and 2-3 zero-collateral knobs — rank 2-5 against deficits of 12-31.  A handle
is the narrowest possible carrier: it has no siblings to cancel against.  The 440 carriers at the
identical signature (3,16,4,12,2,2,16,39017) are handles whose carrier atom sits in the same
16-equation shell, which is the generic shape of the whole class.

## FINAL PICTURE — every carrier class in the instance is now priced
  class                     carriers   best score   any rank > deficit?
  witness placement                1      39026      YES (rank 7 > deficit 4)  <-- unique
  cascade pins                    20      39018      no (rank <= 6, deficit 9-23)
  p-quantised solo handles     1,143      39017      no (rank <= 5, deficit >= 12)
  named non-solo handles           4      39017      no
**The witness placement is the ONLY carrier in the entire instance whose realizable knob-image rank
exceeds its balance deficit.  7 failing is the floor across every carrier class.**

## Step 13 — TWO-SELECTOR RE-PRICING (the single-selector concern, tested)
**The concern was legitimate for the cascade-pin and handle rows, which WERE built from a
one-selector closure.  It was NOT legitimate for the headline row:**

    selectors ON in the 39,026 witness: 2  ->  x_24601 (w-group) and x_2081 (B-side)

So the rank=7 / deficit=4 measurement was never a single-selector artifact.  The witness is a
genuine two-selector state; my table's most important row already lives in the larger space.

**The accumulator, in my frame.**  x_11317 = x_11532 + x_14681 exists here too, and BOTH summands
are FREE INPUTS in my orientation — x_14681 is the very input my cascade closer uses to repair
a9193.  So my frame never "forced x_11532 to zero"; it was free all along, and turning on a second
selector adds no freedom my frame did not already have.  At the witness, x_11532 moves atoms
{38822, 38989, 40772, 41935} and x_14681 moves {9193, 39614} — all OUTSIDE the 8-atom region, so
neither is a zero-collateral knob and neither raises the rank.

`two.py` — two-selector closures priced directly (criterion unchanged, rank > deficit):

  frame                    nz  |R|  |S|  deficit  knobs  rank  failing  score   gap
  ONE-SELECTOR              4   15    6      9      2     2      15    39018     7
  same-u adjacent           8   36   12     24      5     5      36    38997    19
  same-u near              10   33   18     15      2     2      33    39000    13
  same-u mid               10   30   15     15      5     5      30    39003    10
  same-u far                8   36   16     20      2     2      36    38997    18
  same-w adjacent           8   34   13     21      2     2      34    38999    19
  cross u/w                 7   27   15     12      9     7      27    39006     5
  cross u/w 2               7   26   14     12      2     2      26    39007    10

`twosweep.py` — the full cascade-pin sweep re-run from 1-, 2- and 3-selector bases:

  base frame          base score  chain  base deficit/rank  pins: best  maxrank  min(deficit-rank)  wins
  1sel   u[0]            39018     20        9 / 2            39018       5            6            0
  2sel   same-u          38997     22       24 / 5            38997       8           18            0
  2sel   same-u far      38997     31       20 / 2            38998       9           14            0
  2sel   same-w          38999     29       21 / 2            38999       6           17            0
  2sel   witness set     39012     36        9 / 5            39012       9            4            0
  3sel   u+u+w           38986     39       27 / 5            38986      11           17            0
  WITNESS (2 selectors + 4 detached vars)     4 / 7            39026       7          -3          **1**

**RANK DOES RISE with more selectors — maxrank 5 -> 8 -> 9 -> 11 — which is exactly F's intuition,
and my one-selector frame was indeed a smaller space.  But the DEFICIT rises faster: each extra
selector drags in its own pin cascade, adding equations faster than it adds zero-collateral knobs.
The gap never closes: minimum 4, and it is 4 in the witness's own selector frame.  Zero wins in
any frame.  All two- and three-selector closures score WORSE than one selector (38,986-39,012).**

### What actually distinguishes the witness
Not the selector count (2, same as several rows above) but the FOUR DETACHED VARIABLES
x_642, x_28730, x_29854, x_31864.  The witness's selector set closed by my engine reaches only
deficit 9 / rank 5 / 39,012; the witness reaches deficit 4 / rank 7 / 39,026 because those four
variables sit off my closer's manifold.  **The floor is not a selector-count artifact.  What buys
the last 14 equations is the off-manifold placement, and that is what the four detached variables
encode.**
