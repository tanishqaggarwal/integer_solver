# Agent I — RESUME

## HEADLINE
**The instance contains a 256-bit ECDLP on its all-atoms-zero branch.**
Machine-checkable certificate, **15/15 steps PASS**:

```
cd /home/user/integer_solver/solve_lab/agentI_work && python3 certify.py
```
Needs only EQUATIONS.txt (rebuilds its own caches). ~12 min. Writes
`certificate_results.json`; last run's log is `certify.log`.

| step | claim | status |
|---|---|---|
|0| 39,033 eqs, 40,885 distinct atoms, every atom degree <= 2 | PASS |
|1| **3,707 atoms are `X - p*H` and every handle H occurs in NO other atom** => each is exactly `X == 0 (mod p)` with a free quotient => the mod-p abstraction is EXACT | PASS |
|2a| Z-propagation from EMPTY forces 5,624 vars, 0 conflicts; the ONLY large pinned constants anywhere are `p = 2^256-2^32-977` (220 vars) and `K` (9 vars) | PASS |
|2b| mod-p propagation with **effective-support** reduction determines 28,701 vars from 1,156 boolean decisions in ~1 s, leaving EXACTLY 3 violated atoms | PASS |
|2c| those 3 are `X_k - m_k*X35389`, three different m_k, each X_k an independent multiple of X6671 => rank 2 => **X35389 = X6671 = 0 forced** | PASS |
|3| **EXACT symbolic identities** (sympy back-substitution of the instance's own atoms, not numeric agreement): `X35389 = (x2-x1)^2(x3+x1+x2+K) - (y2-y1)^2`, `X6671 = (y3+y1)(x2-x1) - (y2-y1)(x1-x3)`, x1=X12186 y1=X16742 x2=X14853 y2=X24908 x3=X22162 y3=X30213 K=X24453 | PASS |
|4| `u = x + K/3` makes these the STANDARD Weierstrass addition law; the 512 conditional-pin constants become 256 points on `y^2 = x^3 + b` — 219/256 directly, 256/256 after swap/negate; target T on it | PASS |
|5| order **N is a 256-bit PRIME**; the 256 table points form ONE doubling ladder of length exactly 256 | PASS |
|5b| embedding degree > 50, not anomalous (N != p) | PASS |
|6a| every propagation branch point has root set exactly {0,1} — no continuous knob | PASS |
|6b| propagation confluent => A and B are a function of the boolean decisions ALONE | PASS |
|6c| **MECHANISM TEST**: extra selectors RELEASE the accumulator (it is advice, not derived); setting it to the independently computed EC partial sum CLOSES the rung, 4/4 random triples | PASS |
|6c2| **ESCAPE ATTEMPT**: switch OFF the selector pinning (x2,y2) so it is free advice, set it to `T (-) acc` => A = B = 0 outright — and 3 NEW atoms fire that decompile to the same `m*A' + m'*B'` on a different rung | PASS |
|6d| all 512 conditional-pin handles forced to 0 mod p — one free handle would refute the reduction | PASS |
|6e| releasing a coordinate makes A=B=0 solvable (cubic, roots found) and recreates the identical gadget one rung up | PASS |

### Constants for the record
```
p = 2^256-2^32-977                                        (secp256k1 field prime)
b = 64019533680030876408443198762210829058751700634554282185987325820393598524794
    (a = 0; b is NOT 7 — a sextic twist of secp256k1)
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337  PRIME
t = p+1-N = 432420386565659656852420866390673177327
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
G = (31917591553801470078828036568057743875467637605644620066197178005619323650152,
     83364444556352143115103874010002344754157095926378075484791050960431190202517)
T = ((X22162 + K/3) mod p, X30213)
  X22162 = 36200939269128454586076546451607958467047992891178506183612554289882454126226
  X30213 = 44859544763832475231923253825569092119321525945631045653619508440821028887
```
Solving on this branch = finding k with `k*G = T`, k given by the 256 selector bits.
Honest cost: Pollard rho ~ sqrt(pi*N/4) ~ **2^127 group operations**. No MOV/FR
(embedding degree > 50), not anomalous, order prime so no Pohlig–Hellman.

## THE ONE GAP — do not overclaim
"All atoms zero" is SUFFICIENT for all 39,033 equations but NOT NECESSARY. Every
equation is an integer combination of 3–24 atoms; 1,853 atoms occur in exactly one
equation; the 39,026 witness itself has 9 nonzero atoms. Closing the gap needs:
*no nonzero atom vector in the image of the atom map kills every equation.* Open.
Precise status: **the instance contains a 256-bit ECDLP on the all-atoms-zero
branch, and every construction anyone in this lab has run lives on that branch.**

## Independent corroboration that 39,026 is a CODING optimum (new, from the mod-p side)
With the machinery above I can now produce **complete mod-p solutions** (0 violated
atoms, no extra released variables) by cutting only TWO atoms:
free the K pin `a40368` (`X24453 - K`, **1 equation**) and solve `A = 0` for K, then
free one coordinate and solve `B = 0` for it. Verified working cuts (all reach
`A = B = 0`, 0 conflicts, undetermined = baseline 10,047):
```
{a40368, a29331}  13 equations   {a40368, a10066}  14 equations
{a40368, a10067}  15 equations   {a40368, a26748}  15 equations
{a40368, a10065}  15 equations   {a40368, a29334}  15   {a40368, a29333} 15
```
`mincut.py` enumerated every atom in the derivation chain of x1,y1,x2,y2,x3,y3:
**no second cut costs fewer than 12 equations**, so the cheapest mod-p cut is 13
equations => 39,020. `cutlocal.py` (421 atoms on the defect path, singles then
pairs/triples) reports the same floor. The 39,026 witness beats every one of these
because its 7 residual atoms span 12 equations of which **5 cancel** — it is a
cancellation optimum, not a smaller defect. This reproduces prior sessions'
"39,026 is optimal for its residual" from a completely different direction.

`polyroot.py` + the cubic solve also confirm: released coordinates admit exact
solutions of A=B=0 (two roots, one degenerate x1=x2,y1=y2 which makes the gadget
vacuous) — the escape is real but always relocates the check.

## Exhaustive check of the easy part of the ECDLP
Meet-in-the-middle over subset sums of the 256 ladder points: **no k of Hamming
weight <= 4 gives k*G = T** (32,897-point table, 11 s). Weight <= 6 is reachable
(~2.8M table) if anyone wants it; the full problem is 2^127.

## Score status
- Baseline re-verified with `solve_lab/checker.py`: **39,026/39,033**, failing
  [12231,12270,12350,14584,18673,22044,29125]. CONFIRMED.
- My own mod-p re-solve places the residual in 3 atoms spanning 22 equations => 39,011,
  strictly worse than the witness's 7-atoms/12-equations/7-fail placement. Nothing of
  mine beats 39,026; no new best assignment was written.

## Tools (all mine, agentI_work/)
```
parse.py poly.py dag.py   # independent parse -> atoms.pkl / polys.pkl (40,885 atoms)
model.py <assign.json>    # exact scorer, 0.1 s, reproduces checker.py exactly
prop.py                   # exact Z propagation from empty
fp.py boolscore.py        # mod-p propagation, EFFECTIVE-SUPPORT reduction (~1 s/run)
fprun2/3.py loop2.py      # recorded runs, no-good loop over the choice of free inputs
mech.py                   # standalone mechanism test (ladder verification)
cutscan.py                # disable an atom set, re-solve, see if the rest closes
polyroot.py               # F_p root finding (used to solve the A=B=0 cubic)
certify.py                # THE CERTIFICATE
```

## Single highest-value next experiment
The only way left to move the number is a residual placement with more CANCELLATION,
not a smaller one. Concretely: take a complete mod-p state from a 2-atom cut
(`mincut.py` / `cut2_state.pkl`) and, instead of paying the cut's equations outright,
search for a compensating atom vector inside the SAME equations — i.e. solve
`min |{e : sum_a c_{e,a} v_a != 0}|` over atom vectors v in the image of the atom map,
seeded by the cut. That is the one formulation under which 39,026's five cancellations
are not special, and it is also the exact question the certificate's open gap asks.

## Reproduce everything
```
cd /home/user/integer_solver/solve_lab/agentI_work
python3 certify.py                # the certificate, 15/15 PASS, ~12 min
python3 mech.py 11 12             # standalone ladder mechanism test
python3 mincut.py                 # cheapest second cut enumeration
python3 boolscore.py wit          # 1 s mod-p re-solve, prints the 3 conflicts
```
