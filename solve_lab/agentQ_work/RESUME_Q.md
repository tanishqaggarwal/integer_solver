# RESUME_Q — agent Q.  Angle: algebra on the REDUCED problem.

## 0. SCORE
- Baseline **39,026 / 39,033** re-verified by me with `solve_lab/checker.py`
  (`best/new_instance_partial_39026.json`, failing `[12231,12270,12350,14584,18673,22044,29125]`).
- **I did not beat it.**  No infeasibility is claimed — the opposite: see §3.

## 1. THE REDUCTION IS COMPLETE — the whole instance is ONE discrete log
Everything below is measured in `agentQ_work/`, reproducible from cold, and uses only the
integers in the file plus agent F's decode of the stage law.

**(a) The stage law has a conserved cubic.**  F's law is
`l=(b_y-a_y)/(b_x-a_x); out_x=l^2-a_x-b_x-K; out_y=l*(a_x-out_x)-a_y  (mod p)`.
Put `X = x + K/3` (so `3c = K` cancels the offset).  In X the law is the plain chord
construction, whose conserved set is a cubic `Y^2 = X^3 + aX + b`.  Fitting a,b from **two**
leaf pins and testing the other 251 (`inv1.py`):
  **a = 0**,  b = 64019533680030876408443198762210829058751700634554282185987325820393598524794,
  **253/253 leaf pins lie on it** (the other orientation gives 2/253 — the orientation is
  determined, not chosen).  Discriminant `4a^3+27b^2 != 0`, so the cubic is non-singular.

**(b) The law IS an abelian group operation** (`inv2.py`, measured not assumed):
  chordK == the shifted chord law on 198/198 random leaf pairs;
  **associative 297/297, commutative 297/297** on random leaf triples.
  => the 96-stage tree collapses: *the fold of a leaf subset is its GROUP SUM.*  Tree topology,
  stage roles, mux wiring and the 56 undecoded slot pairs are all **irrelevant** to the answer.

**(c) Group order** (`order.py`, computed here by Cornacchia on `4p = L^2+27M^2` for the a=0
case, then *verified* by exact scalar multiplication killing 5 independent points):
  **N = 115792089237316195423570985008687907852837564279074904382605163141518161494337**
  N is a 256-bit **prime** (no factor < 2*10^5, strong-Fermat), `N != p` (not anomalous),
  `p^k != 1 mod N` for k <= 24 (no small embedding degree).  => no Pohlig–Hellman, no transfer.

**(d) The 256 leaves are ONE doubling ladder** (`ladder.py`).  249/253 decoded leaves have their
double also a leaf; the doubling graph is 4 chains of 124/79/41/9 which link tail->head by
exactly one missing doubling each.  Result, **verified exactly**: `L_i == 2^i * G` for
i = 0..255, G = leaf of selector var **x2779**, three exponents (41, 51, 176) being the three
pins F extracted only one constant for.  `ladder.json`.

**(e) The root target is on the same cubic** and is not a leaf.  The two unconditional pins are
the only atoms in the file with a >=60-digit constant of their shape:
  `((x24468 - C1) - x32989)` and `((8863713*(x18956 - C2)) - x14257)`,
  C1 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
  C2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
  T = (C1 mod p, C2 mod p) satisfies the cubic in the (X = x+K/3, Y) chart; (C2,C1) does not.

## 2. THEREFORE
> **EQUATIONS.txt is satisfiable iff `k*G = T` for some k in [1, 2^256-1], and the leaf ON-set of
> a solution is exactly the binary expansion of k.**

## 2b. THE LADDER IS CONFIRMED AGAINST THE RAW FILE (`check3.py`)
The three exponents 41, 51, 176 were *inferred* (F extracted only one constant for their pins).
Predicting their coordinates from `2^i G` and searching the file's own pin constants:
all three predicted x-coordinates **match a literal constant present in the instance**
(exp 41 -> pin selector x33434 wire x28109; exp 51 -> x18184/x32297; exp 176 -> x22579/x24773).
So the ladder `L_i = 2^i G, i = 0..255` is verified against raw instance data, not just inferred.
The deliverable's ON-set {2081, 24601} = ladder exponents **{72, 235}**, i.e. it evaluates
`k = 2^235 + 2^72`, and `kG != T` — which is exactly why it must break 7 atoms.

## 3. FEASIBILITY (positive result, replaces every "coding optimum" reading)
Subset sums of `{2^i G}_{i=0..255}` realise `kG` for **every** k in [1, 2^256-1].  N < 2^256, so
every group element is hit — **including T**.  A satisfying assignment therefore EXISTS.  The
instance is feasible; what is missing is only the *index* k.

## 4. WHAT THIS RETIRES
- **F's "highest-value next experiment" (deep meet-in-the-middle down the 78-side to a support-<=24
  node, enumerate under 2^24) cannot work.**  MITM on a tree only helps if the fold mixes; here the
  fold is a group sum, so *any* MITM is the generic square-root attack: **~2^128**, not 2^24.
  Nothing about the tree shape (178|78 split, the 88- and 50-support stages, the 66 stages in the
  2..24 window) changes that.  Finishing the 56-stage decode buys nothing for the search.
- Likewise the "reachable space is 2^256-1 subsets" count is right but the space is a *group*, so
  it is exactly one coset problem, not an exponential combinatorial search.

## 5. THE DELIVERABLE, RE-READ (val4.py)
Evaluating every atom at `new_instance_partial_39026.json` **as given** (no forward pass):
exactly **7 nonzero atoms, 7 failing equations** — and they are **not** the target pins.
The deliverable *satisfies* the root pins (x13682 = C1, x37892 = C2 mod p) and pays for it by
breaking 7 internal gate/definition atoms.  Its leaf ON-set is `{2081, 24601}` (2 leaves); the
group sum of those two leaves is NOT T, which is exactly why the 7 atoms must absorb the gap.
`mult.py`: every atom sits in 8–20 equations (only one atom in the file sits in exactly 1), so the
7 atoms are chosen so that their contributions **cancel** in all but 7 rows: the 7-atom set touches
12 equation rows and zeroes 5 of them.  A pure forward evaluation of the same free inputs leaves 4
nonzero residual atoms touching 13 rows => 39,020.  `rows4.py` prints both restricted row systems.

## 5b. WHAT THE DELIVERABLE ACTUALLY DOES (measured, mod p, over all 38,748 wires)
Counting wires whose value mod p equals each candidate, in `new_instance_partial_39026.json`:
`L_72 = 2^72 G` x-coordinate on **92** wires, `L_235` on **5**, the group sum `L_72+L_235` on **0**,
the target x-coordinate C1 on **4** (the top).  So the deliverable does **not** fold at all: it
passes a *single* leaf (2^72 G) up the whole tree as a chain of one-live-input pass-throughs, cuts
the second leaf off after 5 wires, and then **overwrites the value with the target near the root**,
paying 7 broken atoms for the overwrite.  That is the entire content of the 39,026 partial.

## 6. WHAT I RAN
| file | result |
|---|---|
| `inv1.py` | cubic invariant fitted; 253/253 leaves on it |
| `inv2.py` | associativity/commutativity 297/297; chordK == group add 198/198 |
| `order.py` | N (256-bit prime), verified by scalar mult |
| `qstruct.py` | 4 doubling chains 124/79/41/9 |
| `ladder.py` | full ladder `L_i = 2^i G`, verified; missing exponents 41,51,176 |
| `val4.py` | deliverable = 7 nonzero atoms / 7 failing; ON-set {2081,24601} |
| `mult.py`,`rows4.py` | atom→equation multiplicities; restricted row forms |
| `dlp_bsgs.py` | BSGS for small k (see `dlp_bsgs.log`) |
| `fastg.py` | Jacobian + gmpy2 group arithmetic (gmpy2 installed by me) |

## 7. NEXT EXPERIMENTS, IN ORDER
1. **Lottery tickets on k** (cheap, each would fully solve): small k (BSGS, running);
   low-Hamming-weight k (meet-in-the-middle over 3-subsets of the 256 ladder points, 2.7M each
   side, covers weight <= 6); k = a*2^s with small a.
2. If those miss, the honest statement is: the remaining work is a **generic 256-bit discrete
   log, ~2^128 group operations**, and no amount of circuit decoding reduces it.
3. Score-side (independent of the DLP): the minimum failing count is
   `min over achievable atom-value vectors r of #{rows with <k,r> != 0}`.  My model adds one
   degree of freedom nobody had: the defect `delta = kG - T` can be made **any** group element by
   choosing k.  But every extra row you try to cancel imposes a *linear* condition on (X,Y) of
   the fold point, i.e. a line meeting the cubic in <=3 points — reaching one of them is again a
   DLP.  So the extra freedom does not lower 7 without solving the DLP.  Measured, not assumed:
   `rows4.py` gives the restricted forms; the deliverable-7 set has a row `1*[3]` (single atom)
   which is what forces its 7th failure.

## 8. STANDING RULES OBSERVED
- No claim of the form "nothing can move X" is made here.
- Everything above is at the level of integers/congruences and the polynomials in the file; I
  neither used nor needed any named-object framing, and I did no generator forensics.
- Verification used `checker.py` on the 39,026 file only (it parses fine); nothing I produced
  needed `verifyE.py`.

## 9. LOTTERY-TICKET SWEEPS ON k  (each would have fully solved the instance)
All use `fastg.py` (Jacobian + gmpy2 + Montgomery batch inversion, ~40k group ops/s/core).

| sweep | family of k covered | result |
|---|---|---|
| `dlp_bsgs.py` | k < 2^44 and N-k < 2^44 | **none** (275 s) |
| `lowwt.py` | Hamming weight(k) <= 6 | **none** (127 s) |
| `wt7.py` | Hamming weight(k) <= 7 | **stopped** at ~2% (CPU contention, load 19 on 4 cores); re-runnable, ~3 h alone |
| `window.py` | all ON-bits inside a 34-bit window (k = a*2^s, a < 2^34) | **none** (2865 s) |
| `smallmul.py` | m*T on the ladder for m <= 10^7, i.e. k = 2^i/m mod N | **none** (249 s) |
| `lam.py` | (1) k = +-lambda^j * 2^i ; (2) k = a + b*lambda, \|a\|,\|b\| < 2^21 | **none, both** (258 s) |

`lam.py` also **confirms the endomorphism**: with beta a cube root of 1 mod p, phi(X,Y) = (beta X, Y)
equals multiplication by a cube root lambda of 1 mod N.  It gives at best a sqrt(3) speedup, so it
does not change the 2^128 figure.

## 10. HONEST BOTTOM LINE
The instance is **feasible** and the remaining work is **exactly one 256-bit discrete logarithm in a
prime-order group with no exploitable structure** (prime order, non-anomalous, no small embedding
degree, only the sqrt(3) CM endomorphism).  Generic cost ~2^128 group operations.  Every structured
short-cut I could construct has been tried and missed.  Nothing about the circuit — the 96 stages,
the 56 undecoded slot pairs, the leaf-support profile, the mux quadrants — reduces that number,
because the fold is a group homomorphism from the selector vector.

The 39,026 deliverable is therefore best understood not as a near-miss on a combinatorial search but
as *the cheapest known way to fake the root value while breaking as few equation rows as possible*.
Improving it is a pure coding problem over the atom incidence matrix (agent A's formulation) and is
independent of the discrete log.

## 11. THE CAVEAT IS CLOSED — AND THE LEAF DECODE IS NOW EXACT (`qpins.py`, `qladder2.py`, `qstages.py`)
Everything in this section is re-derived **directly from `EQUATIONS.txt`**, with no input from any
other agent's directory.  `qextract.py` (a re-run of the shared `extract_atoms.py`) rebuilt the atom
database locally: 32,006 distinct gate atoms, all of the form `+ - *` over integers (no division).

**(a) All 256 leaves, no inference.**  Pin atoms have the shape `(x_g)*((x_w)-(BIGCONST))`.  Scanning
for them gives **256 selectors, each with exactly 2 pins** — the three that earlier had one constant
were an extraction shortfall, not missing data.  With the correct shift (`X = c + K/3 mod p`, K =
`x_24453`, a literal wire) **256/256 leaf points lie on the cubic**.  Doubling closes into a *single*
chain of length 256, and **256/256 satisfy `L_i = 2^i G`**, G = selector `x2779`.  This supersedes
section 2b: nothing is inferred any more.  (`qleaf.json`, `qladder.json`)

**(b) The stage law, at every stage.**  The stage gadget is the division-free chord law:
```
dx = ua-ub    dy = ya-yb
R1 = S*dx^2 - dy^2      S = u3+ua+ub+K        <=>  lambda^2 = u3+ua+ub+K
R2 = A*dx  - B*dy       A = y3+yb, B = ub-u3  <=>  y3+yb = lambda*(ub-u3)
```
(the `+K` is exactly `3*(K/3)`, i.e. the raw-coordinate form of `X3 = lambda^2 - X1 - X2`.)
Searching the atom DAG for this shape finds **383 stage gadgets**: 89 leaf-adjacent, 78 mixed,
216 internal.  Each was tested by **Schwartz-Zippel on random curve points** — random P_a, P_b on the
cubic, (u3,y3) set from the group law, then the *actual* sub-DAG from `EQUATIONS.txt` evaluated:

> **383 / 383 verified, including 89 / 89 leaf-adjacent — and all 383 with orientation (+1,+1),
> i.e. every stage computes the plain sum `P_a + P_b`, with no sign flips.**

None of the 1,532 stage core wires is multi-defined, so the test used the real gate relations.
**The caveat in section 3 is closed**, and closed more broadly than it was stated.

**(c) The gadget census has exactly the shape of a combination tree over 256 leaves.**  Counting how
many of each gadget's four input coordinate wires are hard-wired to 0:

| kind | count | hard-zero inputs |
|---|---|---|
| leaf-adjacent | 89 | 0 — combines two leaves |
| mixed | 78 | 2 — one leaf plus a dummy: a pass-through |
| internal | 191 | 0 — live |
| internal | 25 | 4 — dead |

89 leaf pairs consume 178 leaves; the remaining **78** leaves are exactly the 78 pass-throughs.
That is a binary combination tree over all 256 leaves, and it corroborates `fold = group sum`.

## 12. WHAT IS *STILL* UNVERIFIED (the gap moved, it did not vanish)
Gadget outputs do not feed the next gadget directly — they pass through a **selector/mux layer**
(`u3` wires are consumed by `V*V`, `V+V`, `V-V` gates gated on selector bits).  I verified the
*law* each gadget enforces as a function of its four input coordinate wires; I did **not** verify
that the selector logic can realise an arbitrary subset of leaves.  So the existence result
"a satisfying assignment exists" now rests on the **routing layer**, not on the stage law.
This is a strictly smaller and different gap than the one I flagged, but it is real and it is
load-bearing.  Do not report the existence result as unconditional.

## 13. ON THE CROSS-MODEL TENSION (RETRACTED BY ITS AUTHOR - no longer live)
A residual-side exhaustion result briefly appeared to contradict the existence claim.  Its author
retracted it: the closure was **base-local**, not global (its base configuration had two selectors
already on, so a single-selector configuration was three flips away, and configurations generated
outside the BFS's reach landed outside the closed image).  Nothing in it forbids a satisfying
assignment.  I record the retraction and take no credit for it; section 11 was measured on its own
merits and does not depend on this.  For the record, the measurement I ran anyway:

My model predicts the root value is essentially unconstrained: for 300 random **weight-128**
configurations, the fold equalled `k*G`, was on the curve, and gave **300 distinct values**
(300/300).  Under the fold model the reachable set of root values is all of Z/N, so no genuine
*global* closure at a few dozen tuples can exist.  Section 12 still stands on its own: the remaining
unverified link in my chain is the **routing layer**, not the stage law.

## 14. THE ROUTING TEST — RUN NON-CIRCULARLY, AND IT DOES NOT CLOSE (`qsolve.py`, `qrun2.py`, `qdegen2.py`)
`qsolve.py` parses **every term of every equation** (47,198 distinct terms) and does unit
propagation mod p: any term with exactly one unknown variable is solved for it (linear, or
quadratic with a repeated root).  Nothing is assigned by hand, so nothing is presupposed.

**(a) Selectors only — the non-circular run (`qrun2.py`).**  Set the 256 selector bits and nothing
else; let the leaf coordinates be *solved* from the pin atoms.

| weight | ON-leaf X solved | OFF-leaf X solved | OFF-leaf Y forced to 0 | gadget outputs | root |
|---|---|---|---|---|---|
| 1 | **0/1** | 0/255 | 219/255 | 25/383 | no |
| 2 | **0/2** | 0/254 | 219/254 | 25/383 | no |
| 3 | **0/3** | 0/253 | 217/253 | 25/383 | no |
| 5 | **0/5** | 0/251 | 215/251 | 25/383 | no |
| 7 | **0/7** | 0/249 | 214/249 | 25/383 | no |
| 128 | **0/128** | 0/128 | 111/128 | 13/383 | no |

At weight 1: 520/8,583 free inputs solved, 3,236/38,748 wires known, **0 contradictions**.
Turning a selector ON does **not** put that leaf's coordinate on any wire.  This reproduces agent
T's 0-of-256 arrival count independently, from a different parse, and it **contradicts** the report
that liveness is fully determined by the selectors with a configuration space of exactly 2^256.
**In my frame too: routing is a constraint, not a propagation.**  "Set the selectors and evaluate"
is not a well-posed test of this instance.

**(b) OFF leaves are not the identity as points.**  With a selector OFF the leaf's y-wire is forced
to **0** (measured at ~86% of OFF leaves; the rest stay unknown).  The group has prime odd order, so
there is no 2-torsion and `(w,0)` is not a curve point.  Identity behaviour therefore cannot come
from the leaf value — it must come from the **mux coefficients**.  That layer is the crux and I have
not verified it.

**(c) The degenerate branch is vacuous, not doubling (`qdegen2.py`).**  Feed a gadget two EQUAL live
points.  Then dx = dy = 0, so R1 = S*dx^2 - dy^2 = 0 and R2 = A*dx - B*dy = 0 **whatever the output
is**: with the output set to a random wrong value the residual still vanishes, **383/383**.  So the
circuit does not implement point doubling; where two coinciding points meet at a gadget the output
is unconstrained.  The fold picture needs no two equal points ever to meet, which is not guaranteed.

## 15. RETRACTION OF THE STANDING OF SECTION 9
`dlp_bsgs.py`, `lowwt.py`, `wt7.py`, `window.py`, `smallmul.py`, `lam.py` all computed the fold
**inside the group model** and never checked that the circuit agrees at those configurations.
Section 14(a) shows the circuit-side check does not close, and the low-weight regime (1,2,3,5,7) is
exactly where those sweeps live.  **Their clean-miss verdicts are therefore evidence about the group
model, not about the instance**, and I withdraw them as instance-level evidence until the mux layer
is verified.  The searches themselves are correct and re-runnable; it is their standing that changes.

## 16. STATUS OF THE EXISTENCE RESULT
Still **conditional**, and now on a sharper thing than in section 12: it holds if and only if the
**mux-coefficient layer** makes an OFF leaf act as the identity and routes each gadget's output to
the next gadget's input.  Sections 11(b) and 11(c) are unaffected and remain measured: all 383
gadgets enforce plain `P_a + P_b` for distinct inputs, and the gadget census is a combination tree
over 256 leaves.  What is *not* established is that the selector bits pick out a subset at all.

## 17. THE MUX LAYER, SOLVED SYMBOLICALLY AT ONE SLOT (`qmux.py`, `qquad.py`)
Not propagated through — read off `EQUATIONS.txt` as atoms and then checked numerically.
Slot: inputs leaf `2^0` (selector `x_2779`) and leaf `2^164` (selector `x_34715`), chord output
`(x_22294, x_33676)`.  Verbatim atoms:
```
x_2779*(x_2779-1)            x_34715*x_34715-x_34715      <- both selectors boolean-pinned
x_3565 = a                   x_31966 = b
x_24678 = 1-b                x_24849 = 1-a
cA = x_13201 = x_3565*x_24678   = a(1-b)
cB = x_33391 = x_31966*x_24849  = b(1-a)
cC = x_4639  = x_31966*x_3565   = a*b
Xout = x_20820 = cA*x_22231 + cB*x_11321 + cC*x_22294
Yout = x_18440 = cA*x_27051 + cB*x_37031 + cC*x_33676
live_out = x_11830 - x_1609 = (a+b) - ab = a OR b
```
Evaluated on the real leaf constants (`qquad.py`), all four quadrants:

| (a,b) | cA | cB | cC | slot carries | matches |
|---|---|---|---|---|---|
| (0,0) | 0 | 0 | 0 | identity (0,0) | **yes** |
| (1,0) | 1 | 0 | 0 | leaf 2^0 | **yes** |
| (0,1) | 0 | 1 | 0 | leaf 2^164 | **yes** |
| (1,1) | 0 | 0 | 1 | **sum 2^0 + 2^164** | **yes** |

**The mux does implement identity, pass-through and sum.**  This is L's mutually-exclusive-quadrant
claim, confirmed in my own frame rather than taken.  It also resolves the §14(b) worry: the identity
value is `(0,0)`, not a curve point, but it is only ever *passed through* — it can never enter a
chord, because `cC = ab = 0` whenever a child is dead.

**Generality.** An automated association-free structural match of this exact shape succeeds at
**188 / 383** slots (19 with both selectors boolean-pinned, 169 with internal live bits).  The other
195 carry the same `c*u3` product but my matcher did not confirm their summation tree; they are
*consistent with* the law, **not confirmed by this test**.  One slot is done completely; the rest is
a matching problem, not a semantic one.

**Why §14(a) stalled — now explained.**  A leaf pin is not `sel*(w-C)`; it is `sel*(w-C) - z` for a
further wire z.  The coordinate is forced onto the wire only once z is separately forced to 0, so
unit propagation from the selectors alone can never place a leaf coordinate.  The routing is
determined, but by a *simultaneous* system, not by propagation — exactly as agent T said.

## 18. WHERE THE EXISTENCE RESULT NOW STANDS
It closes **if** (a) the quadrant law of §17 holds at all 383 slots (confirmed at 188), and (b) no
two *equal* points ever meet at a live slot, where §14(c) showed the chord residual is vacuous.
For (b) there is a clean criterion: children of a slot are sums over **disjoint** leaf subsets, so
they coincide only if `sum_{S1} 2^i - sum_{S2} 2^i = ±N`.  Both sums are `< 2^256 < 2N`, so no other
collision is possible.  That is a checkable condition on the particular k, not a generic hazard.
§15 stays in force: the §9 sweeps regain instance-level standing only when (a) is closed at 383/383.

## 19. (d) CLOSED — THE QUADRANT LAW HOLDS AT 383 / 383 (`qmux2.py`, `qtree.py`, `qlivetree.py`)
The 195 "unmatched" slots of §17 were a **labelling** artefact, not a different law: `qstages.py`
picked `(u3,y3)` as `sorted(free)[0],[1]`, which is the X/Y order at some slots and the reverse at
others.  Letting the matcher try both orders, and requiring **both coordinate muxes to use the
identical coefficient wires** `cA,cB,cC` (which rules out accidental structural matches):

> **383 / 383 slots confirmed** — 40 with both live bits boolean-pinned, 343 with internal live bits.
> Zero unmatched.

**The liveness composition is a single tree (`qlivetree.py`).**
* every one of the 383 slots emits `OR(s1,s2)` of its own two live bits — **383/383**;
* those ORs give **382 parent<-child edges** among 383 slots, i.e. exactly a tree;
* **exactly one slot has no parent**, and **all 383 are reachable from it**;
* **all 256/256 leaf selectors appear under that single root**;
* the 766 live-bit slots decompose as **256 leaf selectors + 382 child ORs + 128 hard zeros**
  (the 128 zeros are the dead dummy branches of the pass-through slots — §17's `cB = cC = 0` case).

So the slots are one tree over all 256 leaves, each node computing identity / pass-through-left /
pass-through-right / sum according to `cA = s1(1-s2)`, `cB = s2(1-s1)`, `cC = s1 s2`.
**On the liveness side the fold picture is now measured end to end, from the 256 boolean leaf
selectors to a single root.**

## 20. THE ONE THREAD STILL LOOSE — COORDINATE HAND-OFF
A slot's mux output wires are **not literally** the next slot's input wires: `qtree.py` finds
0 / 383 slot outputs feeding another slot directly, and the root pin `x_24468` is not a slot output
but `x_24468 = x_13682 + 12354891 * x_34243`.  There is an additional additive/aliasing layer
between a slot's mux output and its parent's input, and between the top slot and the root pin.
The liveness tree is isomorphic to the point tree and the quadrant law ties coordinates to live bits
at every slot, so the coordinate composition almost certainly follows the same tree — but
**I verified the tree on the liveness side, not on the coordinate side**, and I am not going to
claim the coordinate hand-off on the strength of an isomorphism I did not measure.
Existence therefore remains conditional on that one layer; §15 stays in force.

## 21. (e) CHARACTERISED — IT IS AN AFFINE ALIAS, BUT THE SLACK IS NOT PINNED (`qalias.py`)
What sits between a slot's mux output M and its parent's input P, read verbatim off the instance:
```
x_17675 - x_20820 - x_36780             P = M + Q        Q = x_36780 = x_4116 * x_22163
6910381*(x_15439 - x_18440) - x_11630   k*(P - M) = Q    Q = x_11630 = x_1962 * x_10858
x_24468 - x_13682 - 12354891*x_34243    ROOT PIN = M + k*Q,  Q = x_34243 = x_16153 * x_14393
```
Across all 766 mux outputs: **573 alias to a parent slot input, 2 alias to the ROOT PIN**, 191 use a
shape I did not chase.  Forms: `P = M + Q` 192, `k*(P-M) = Q` 192, `P = M + k*Q` 191.
**All 575 slack wires are products of two wires.**  So layer (e) *is* the affine alias
`parent_input = mux_output + (multiple of) a product`, and the root pin is the top slot's mux output
under exactly that alias — the coordinate composition follows the same tree measured in §19.

**But the slack is not forced to zero, as far as I can measure.**  Of the 575 slack products,
523 have both factors used elsewhere and 52 have a factor occurring in exactly ONE term (wholly
unconstrained).  The shared factors — `x_4116` (66 terms), `x_16153`, `x_1962`, `x_12682`,
`x_19049`, `x_15616` — carry **no unary pin at all**: no boolean constraint, no zero pin.  I could
not exhibit anything forcing `Q = 0`.  Until something does, the parent input is the mux output plus
an unpinned amount, and **the coordinate hand-off is not determined by the tree**.
So (e) does **not** close, and **§15 stays in force**.

**Bearing on agent K's null result.**  Two candidate explanations, and I measured only the first:
1. **Measured (my §5b):** in the 39,026 deliverable the group sum appears on **0** wires — because
   that assignment never folds at all (one leaf propagates 92 wires, the other is cut after 5).
   At any configuration that does not fold, a search for the composition must come up empty
   whatever the aliasing does.
2. **Structural, not measured:** if the alias slack is nonzero, the composition sits on the mux
   output wire but *not* on the parent's input wire, so a literal search would find it on at most
   one wire even in an assignment that does fold.
K's null therefore does not by itself show the circuit fails to force compositions; but neither does
my work show that it does.  I would not let a barrier withdrawal rest on the null alone.

## 22. RULING ON SECTION 15 — PARTIAL RESTORATION, WITH A NEW AND DIFFERENT GATE (`qmult.py`)
L's result closes the hand-off mod p: the six shared factors are the constant p (220 such wires in
the instance), so every slack term is `p * (free variable)` and `parent_input = mux_out (mod p)`
unconditionally.  My §21 could not see this because I looked for something *forcing* the factors to
zero; nothing needs to force a constant.  The 523/52 split I reported was one population, not two.

**Does the mod-p qualifier hurt the six searches?  No — and here is why.**  `dlp_bsgs.py`,
`lowwt.py`, `wt7.py`, `window.py`, `smallmul.py`, `lam.py` are all *negative* results of the form
"no k in family F has kG = T".  The direction they need is

>  a satisfying assignment with ON-set S exists  ==>  k = sum_{i in S} 2^i satisfies kG = T

and every step of that implication is a point identity **mod p** — the leaves, the chord gadgets,
the quadrant muxes, the tree, and now the hand-off.  The 927 `c > 1` conditions over Z are about
**constructing** an integer lift, i.e. the *converse* (existence) direction.  So the Z gap sits on
the existence claim, not on the negatives.  **On that count mod p is exactly the modulus they
needed**, and my §15 condition is met.

**But measuring it exposed a different assumption I had never checked, and it does not hold on rank
grounds.**  I had been assuming that satisfying the bundled equations forces the individual atoms to
zero.  Counting occurrences **without** deduplication (the earlier "every atom in 1 equation" figure
was an artefact of `gates.jsonl` being deduplicated):

* **47,198 distinct atom terms across 39,033 equations**, mean **11.5 atoms per equation**;
* **82.7%** of atoms occur in **>= 2** equations (mode 11, max 19); 8,166 occur in exactly one;
* stage-core atoms: 275 in one equation, the rest spread over 4-19.

The incidence matrix therefore has **more columns (47,198) than rows (39,033)**, so its null space
has dimension **>= 8,165**.  **The bundling does not force the atoms to zero on rank grounds.**
That does not prove compensation is achievable — atoms are constrained functions of wires, not free
coordinates — but it does mean the implication above is not established.

**Ruling.**  The six programs move from *group-model only* to **instance-level conditional on
atom-forcing** — strictly stronger than where §15 left them, strictly weaker than unconditional.
I am not restoring them fully.  The gate that remains is new, smaller, and precisely stated:

> Does the null space of the atom-incidence matrix contain a vector **realisable by an actual wire
> assignment**?  Rank counting says the null space is non-empty; realisability is the open part.

That is agent A's object, not mine.  It also bears on the lab's scoring frame, which treats the atom
as the unit of failure: the 39,026 partial has 7 nonzero atoms and exactly 7 failing equations, yet
the median atom sits in 11 equations — worth checking whether its 7 lie among the 8,166 singletons.

**Follow-up measurement (same script).**  The 39,026 partial's 7 failing equations
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]` contain 20, 8, 24, 20, 3, 2 and 15 atoms.
Only **eq 22044** contains singleton atoms (both of its 2).  In the other six, every atom occurs in
**6-15** equations.  So a single nonzero atom there would generically break many equations, yet only
7 break in total.  Either the nonzero atoms are very few and concentrated, or **compensation between
atoms is already happening in the lab's best assignment**.  That is a concrete, cheap handle on the
atom-forcing question and it belongs with agent A's incidence matrix.

## 23. RECONCILIATION WITH F'S PARSE — MY RANK ARGUMENT WAS AN ARTEFACT, AND §22'S GATE IS CLOSED
Measured against my own `qmult.pkl`:

| | count |
|---|---|
| my distinct terms | **47,198** |
| ...occurring in **>= 2** equations | **39,032**  <- F reports **39,033** atoms |
| ...occurring in **exactly 1** equation | **8,166** |
| my excess over F (47,198 - 39,033) | **8,165** |

**My multi-occurrence terms match F's atom count to within 1, and my 8,166 singletons account for
exactly the 8,165 excess.**  So the entire "nullity >= 8,165" observation of §22 is a restatement of
my own extra granularity: my parser splits sub-expressions that F treats as single atoms, each split
landing in exactly one equation.  It was never a statement about the instance.  I said in §22 that
my finer atoms are constrained functions of wires rather than free coordinates; that caveat turns
out to be the whole story.

With F's `ker(M) = 0` (three independent computations) and T's faithfulness check (exact list
equality between `{e : (Ma)_e != 0}` and `checker.evaluate_all` at 10 points), all-atoms-zero is an
**equivalence** in the 39,033 decomposition: any assignment satisfying all 39,033 equations makes
every atom vanish.  **The atom-forcing gate I opened in §22 is closed.**

## 24. FINAL RULING ON SECTION 15 — FULL RESTORATION
The implication the six negatives need,
`satisfying assignment with ON-set S  ==>  k = sum_{i in S} 2^i satisfies kG = T`,
now has every link measured: all-atoms-zero is forced (§23); the leaves are `2^i G`; 383/383 chord
gadgets compute `P_a + P_b`; 383/383 slots implement identity / pass-through / sum; the slots form
one tree with a single root over all 256 leaf selectors; and the coordinate hand-off follows that
tree **mod p** — which §22 established is exactly the modulus this direction needs, the 927 `c > 1`
conditions over Z sitting on the *converse* (existence) direction.

> **`dlp_bsgs.py`, `lowwt.py`, `window.py`, `smallmul.py`, `lam.py` are restored to full
> instance-level standing.  `wt7.py` is restored at its true coverage: 33.7% of the weight-<=7
> sweep, no hit in that portion — a partial, not a bound.**

The weight-<=6 bound, the k < 2^44 bound, the 34-bit-window bound, the m <= 10^7 bound and the
lambda-basis bound are statements about **the instance**, not about the group model.

## 25. (d) SETTLED — THE ATOM IS NOT THE UNIT OF FAILURE
`ker(M) = 0` forbids *all* equations being satisfied with some atom nonzero.  It does **not** forbid
an atom being nonzero in an equation that still sums to zero.  The 39,026 partial shows exactly
that: its 7 failing equations contain 20, 8, 24, 20, 3, 2 and 15 atoms, and outside eq 22044 every
atom in them occurs in **6-15** equations — so most occurrences of those nonzero atoms **cancel**,
and only 7 equations break.  **Compensation between atoms is already happening in the lab's best
assignment.**  Scoring as though one bad atom costs one equation is therefore not exact, and the gap
runs in the favourable direction: an atom can be wrong in many equations and cost only a few.
