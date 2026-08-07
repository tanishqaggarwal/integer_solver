# Agent D — resume brief (FINAL)

## Verdict
**The instance is a 256-bit ECDLP on a curve isomorphic to secp256k1, with prime group order.**
Derived independently, from the residual downward, without importing anyone's artifacts.
A full solve is therefore not reachable by search; the deliverable stays `39,026`.

## Best verified
* Global best unchanged: `solve_lab/best/new_instance_partial_39026.json` = **39,026** (re-verified by me).
* My own best independent state: `agentD_work/D_39017.json` = **39,017**, checker-verified,
  failing [56,133,2071,8073,13660,15299,16622,17726,19066,22093,25480,28653,31061,32894,34517,34892].
  Residual = only THREE atoms {a688, a1618, a40608} (the two mod-p pins on the output point).

## The reduction, in five lines
1. Forward-eval from free inputs kills all 31,475 gate atoms; only checks can fail.
2. 192 "advice" congruences `K*(u_free - w) - p*handle` — one Gauss–Seidel sweep (`adv3.py`) closes
   them all: **39,013** (`D_adv.json`). Residual is then only `sel * (combo of A,B) == 0 (mod p)`.
3. `A = x_35389 = (x1+x2+x3+K)(x2-x1)^2 - (y2-y1)^2`, `B = x_6671 = (y3+y1)(x2-x1) - (y2-y1)(x1-x3)`,
   verified digit-for-digit. The whole instance holds iff A ≡ B ≡ 0 (mod p).
4. K = x_24453 is the `a2` of a general Weierstrass curve. Fitting a4,a6 → **256/256** gated
   table points and the pinned target lie on it. Depressed form: **y^2 = x^3 + B, A = 0**,
   B = 64019533680030876408443198762210829058751700634554282185987325820393598524794,
   B/7 is a 6th power mod p ⇒ isomorphic to secp256k1; order = n_secp256k1 (prime).
5. The 256 gated points are one **doubling chain** of length 256 rooted at bit x_2779, so the
   assertion is `k·P_0 = T` — a 256-bit discrete log. `ecdlp.json` has B, shift, ladder, T.

## Pipeline (self-contained in agentD_work/, no s9/s10 caches)
```
cd solve_lab/agentD_work
python3 build_cache.py                  # cache/{atoms,polys,gates,topo}.pkl  (~45 s)
python3 adv3.py D_state1.json D_adv.json # 39,002 -> 39,013
python3 ecsolve2.py D_adv.json           # CRT solve of A,B -> 39,009, then handles -> 39,017
python3 gens26b.py && python3 lat26b.py  # proves MAX 5 of 12 at the 39,026 placement
```
Tools: `dlib.py` (exact eval), `engine2.py` (**St with block=** — incremental apply/revert, ~5 ms),
`rad.py` (reverse AD mod 2^61-1), `hsweep.py`, `scanAB.py`, `scanpairs.py`, `table.py`, `banks.py`,
`condpins.py`, `intsolve.py` (column-HNF integer solver), `placement.py`.

## TRAP that will bite anyone searching from the 39,026 witness
Its five deliberately-broken **gate** atoms {22229, 22230, 35758, 35761, 35762} are silently
repaired by any plain forward ripple, which collapses the score to ~39,008. Always construct the
state as `E.St(v, block={22229,22230,35758,35761,35762})`.

## Optimality, re-derived independently
With gates blocked, the cost-free generator lattice at the 39,026 placement is **9 generators**
(x_642, x_1329, x_8731, x_9118, x_9413, x_10903, x_17325, x_29854, x_31864); only a22231 is frozen.
Enumerating all 2^12 subsets of its 12 equations with an exact integer (HNF) solver:
sizes 12..6 all infeasible over Z, size 5 feasible ⇒ **failing = 7 ⇒ 39,026 is exact.**

## Single highest-value next experiment
None inside the search paradigm. The only thing that changes the answer is the discrete log
`k·P_0 = T` on secp256k1 — so the next experiment is a *feasibility* one: verify the ladder root
P_0 against the standard secp256k1 generator under the isomorphism (x,y) -> (x/u^2, y/u^3),
u^6 = B/7, and confirm T is not a small or otherwise special multiple of P_0 (BSGS to 2^40,
a few core-hours) before declaring the instance closed at 39,026.
