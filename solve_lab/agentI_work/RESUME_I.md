# Agent I — RESUME (build-from-scratch / complete-search angle)

## HEADLINE RESULT (new, checkable, changes the whole picture)
**The instance reduces exactly to a 256-bit ECDLP on a prime-order elliptic curve.**
Derived from scratch, not assumed. Chain of evidence, each step reproducible:

1. Propagation from EMPTY over Z (`prop.py`) forces 5,624 vars with zero conflicts.
   The ONLY large constants pinned anywhere are `p = 2^256-2^32-977` (220 vars,
   secp256k1's field prime) and `K = 97553848499418123410591666447050222001188385549510401465815187079080512838891`
   (9 vars). Everything else pinned is 0 or 1.
2. **3,707 atoms have the form `X - p*H` where every handle H occurs in NO other atom**
   (verified). So each is exactly the assertion `X == 0 (mod p)` with a free quotient
   => the mod-p abstraction of this instance is EXACT.
3. Mod-p propagation with *effective-support* reduction (a var whose coefficient dies is
   not an unknown) determines **28,701 / 38,748 vars from 1,156 boolean decisions in ~1 s**
   (`boolscore.py`). Exactly **3 atoms** are then violated:
   `a17810: X2287 - 8272701*X35389`, `a17813: X21889 - 8646263*X35389`,
   `a17816: X25156 - 10159099*X35389`. Rank 2 in (X35389, X6671) => forces
   **X35389 = X6671 = 0**.
4. Tracing those two (independently reproduced):
   `X35389 = (x2-x1)^2*(x3+x1+x2+K) - (y2-y1)^2`,
   `X6671  = (y3+y1)(x2-x1) - (y2-y1)(x1-x3)`,
   with x1=X12186, y1=X16742, x2=X14853, y2=X24908, x3=X22162, y3=X30213. Verified
   digit-for-digit against the propagated values.
5. **The K offset is removable**: substituting u = x + K/3 turns these into the standard
   short-Weierstrass addition law. Under that substitution the 512 conditional-pin
   constants become points on **y^2 = x^3 + b**,
   `b = 64019533680030876408443198762210829058751700634554282185987325820393598524794`
   (219/256 directly, the other 37 after swapping/negating the pin pair) — a sextic
   twist of secp256k1 (same p, a = 0, different b). The target point is on it too.
6. **Group order N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
   — 256 bits and PRIME** (Cornacchia 4p = L^2+27M^2, order verified by killing a point).
7. **185 of the 219 on-curve table points have their double also in the table** => the 256
   selectors index a doubling ladder {2^i G}. The assertion the instance makes is
   `SUM_{i: b_i=1} 2^i G = P_target`, i.e. **k*G = P_target with k the 256 selector bits**.
8. All 512 conditional-pin handles are forced to 0 mod p (checked), and every "door"
   factor (X38100, X22399, X23917, X11360, ...) is a copy of p, i.e. 0 mod p. So there is
   **no mod-p knob anywhere except the selector bits**, and setting them is exactly ECDLP.

=> **No repair of the 39,026 basin, and no search over this instance, reaches 39,033
without solving a 256-bit ECDLP.** 39,026 is a coding optimum, not a near-miss.

## Status / best score
- Baseline re-verified by me: 39,026/39,033 (`best/new_instance_partial_39026.json`).
- My own mod-p re-solve places the residual in 3 atoms spanning 22 equations => 39,011.
  WORSE placement than the witness's (7 atoms, 12 equations, 7 fail). Not written out.

## Commands
```
cd /home/user/integer_solver/solve_lab/agentI_work
python3 parse.py && python3 poly.py && python3 dag.py   # caches (~1 min)
python3 model.py <assign.json>     # exact score, 0.1 s (matches checker.py)
python3 boolscore.py wit           # 1 s mod-p re-solve, prints the 3 conflicts
```

## Next experiment
Test the curve for the only remaining weaknesses (embedding degree / MOV, small k by
BSGS). If none, the ceiling is a placement question: minimise #equations whose core is
nonzero, over reachable residual vectors. My mod-p machine can enumerate distinct
placements cheaply (one per boolean setting) — look for any placement < 7.
