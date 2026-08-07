# Agent I — RESUME

## HEADLINE: the instance contains a 256-bit ECDLP on its all-atoms-zero branch

**Machine-checkable certificate: `python3 solve_lab/agentI_work/certify.py`**
(needs only EQUATIONS.txt; rebuilds its own caches; prints PASS/FAIL per step;
writes `certificate_results.json`; log in `certify.log`).

Steps 0–5b VERIFIED PASS. Step 6 is the adversarial/weight-bearing step.

| step | claim | status |
|---|---|---|
|0| 39,033 eqs, 40,885 distinct atoms, every atom deg<=2 | PASS |
|1| **3,707 atoms are `X - p*H` and every handle H occurs in NO other atom** -> each is exactly `X == 0 (mod p)` with a free quotient, so the mod-p abstraction is EXACT | PASS |
|2a| Z-propagation from EMPTY forces 5,624 vars, 0 conflicts; only large pinned constants anywhere are `p=2^256-2^32-977` (220 vars) and `K` (9 vars) | PASS |
|2b| mod-p propagation with effective-support determines 28,701 vars from 1,156 boolean decisions, leaving EXACTLY 3 violated atoms `X_k - m_k*X35389` | PASS |
|2c| three different m_k, each X_k an independent multiple of X6671 -> rank 2 -> **X35389 = X6671 = 0 forced** | PASS |
|3| **exact symbolic identities** (back-substitution of the instance's own atoms, sympy, not numeric agreement): `X35389 = (x2-x1)^2(x3+x1+x2+K)-(y2-y1)^2`, `X6671 = (y3+y1)(x2-x1)-(y2-y1)(x1-x3)` with x1=X12186 y1=X16742 x2=X14853 y2=X24908 x3=X22162 y3=X30213 K=X24453 | PASS |
|4| `u = x + K/3` turns these into the STANDARD Weierstrass addition law; the 512 conditional-pin constants become 256 points on `y^2 = x^3 + b`, 219/256 directly and 256/256 after swap/negate; target T on it | PASS |
|5| group order **N = 115792089237316195423570985008687907852837564279074904382605163141518161494337 is a 256-bit PRIME**; the 256 table points form ONE doubling ladder of length exactly 256 | PASS |
|5b| embedding degree > 50, not anomalous (N != p) | PASS |
|6a-e| adversarial: no non-boolean branch; propagation confluent; **mechanism test** — circuit's (x1,y1) equals the independently computed EC sum of the switched-on ladder points; all 512 pin handles forced 0 mod p; releasing a coordinate makes A=B=0 solvable but re-creates the identical gadget one rung up | running |

**Constants for the record**
```
p = 2^256-2^32-977                (secp256k1 field prime; curve b is NOT 7)
b = 64019533680030876408443198762210829058751700634554282185987325820393598524794
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337  (prime)
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
G = (31917591553801470078828036568057743875467637605644620066197178005619323650152,
     83364444556352143115103874010002344754157095926378075484791050960431190202517)
T = (x3 + K/3 mod p, y3), x3=X22162=36200939269128454586076546451607958467047992891178506183612554289882454126226
                          y3=X30213=44859544763832475231923253825569092119321525945631045653619508440821028887
```
Solving the instance on this branch = finding k with k*G = T. Honest cost:
Pollard rho ~ sqrt(pi*N/4) ~ **2^127 group operations**. No MOV, no Smart, no smoothness.

## THE ONE GAP (stated plainly — do not overclaim)
"All atoms zero" is SUFFICIENT for all 39,033 equations but NOT NECESSARY: each
equation is an integer combination of 3–24 atoms, 1,853 atoms occur in exactly one
equation, and the 39,026 witness itself has 9 nonzero atoms. Closing the gap needs
the compensation-closure result (no nonzero atom vector in the image kills every
equation) — that is open. So: **the instance contains a 256-bit ECDLP on the
all-atoms-zero branch, and every construction anyone has run lives on that branch.**
That is why the 39,026 floor has never moved.

## Score status
- Baseline re-verified by me with `checker.py`: 39,026/39,033, failing
  [12231,12270,12350,14584,18673,22044,29125]. CONFIRMED.
- My own mod-p re-solve puts the residual in 3 atoms spanning 22 equations => 39,011.
  Strictly worse placement than the witness's (7 atoms / 12 equations / 7 fail).
  Not written out; nothing of mine beats 39,026.

## Tools (all mine, in agentI_work/)
```
parse.py poly.py dag.py     # independent parse -> atoms.pkl / polys.pkl
model.py <assign.json>      # exact scorer, 0.1 s, reproduces checker.py exactly
prop.py                     # exact Z propagation from empty
fp.py boolscore.py          # mod-p propagation with EFFECTIVE-SUPPORT reduction (1 s)
fprun2.py fprun3.py loop2.py# recorded runs + no-good loop over free-input choice
cutscan.py                  # disable an atom set, re-solve, see if the rest closes
certify.py                  # THE CERTIFICATE
```

## Single highest-value next experiment
`cutscan.py`: for every atom in <=6 equations, disable it (allow it nonzero) and
re-run the mod-p solve. Any atom whose removal makes the rest consistent absorbs the
whole defect at a cost equal to its equation count — a hit below 7 beats 39,026 and
would also be the first real evidence about the compensation gap above.
