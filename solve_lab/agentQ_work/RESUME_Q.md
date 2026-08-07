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
