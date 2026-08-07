# Agent C — RESUME.  Best verified of mine: 39,013 (agentC_work/BEST_39013.json).

Setup: rebuild caches `cd solve_lab/s9 && python3 atomize.py && poly.py && gates.py && fwd.py`.
Installed (absent before): z3-solver, python-sat, cvc5, python-flint, ortools, sympy, numpy.

## Established (independent parse; scripts all in solve_lab/agentC_work/)
1. Break the 900 duplicate-equality 2-SCCs -> pure DAG, **8,173 free inputs**, 30,575 gates all with
   output coefficient 1.  Forward from free inputs = 0: 0 broken gates, **6 nonzero checks**, 39,005.
2. Those 6 = 3 conditions: `x_18956 = K1 (mod p)`, `x_24468 = K2 (mod p)`, `OR(s1,s2) = 1`.
   In branch (1,1) the outputs are FREE inputs x_22162, x_30213 -> setting them to K2, K1 makes all
   three conditions EXACTLY zero over Z (DECISIVE.py).  Closing the two activated bits' pins -> 39,013.
3. Residual = the EC point-addition law.  Curve `y^2 = x^3 + a2 x^2 + a4 x + a6` (curve.json),
   j = 0, short form A = 0, **group order = n_secp exactly** => isomorphic to secp256k1.
   Q = (K2 mod p, K1 mod p) is on it.  All 256 free leaf bits carry pinned curve points forming a
   **doubling chain P_i = 2^i G** (verified i = 1 and i = 255); 178 exponents on the s1 side, 78 on s2.
4. Two doors: `P1 + P2 = Q` (ECDLP, dead; instance literals tested as the dlog, no hit) or
   `P1 = P2` (A, B vanish identically, root output free).  **P1 = P2 is EXACTLY IMPOSSIBLE**:
   it forces kA - kB = +-n with disjoint supports and both deterministic carry chains overflow
   (carry2.py).  The 39,026 deliverable FAKES it by detaching leaf-pin handle atoms.
5. Cost model (mincost.py, cluster.py): faking P1 = P2 costs `|E(hx) u E(hy)|` for the overridden
   bit; min over all 256 bits = **7, at bit x_10513** (a8427: 7 eqs, a8429: 5 eqs, union 7) vs 11 for
   the deliverable's x_24601.  Cluster balance: x_10513 gives |E| = 7, n = 3, |E| - n = 4 — the same
   slack as the deliverable's 12/8 cluster in half the equations, so `failing = 7 - 3 + c`
   (6 -> 39,027 if c = 2; 4 -> 39,029 if c = 0).

## SINGLE HIGHEST-VALUE NEXT EXPERIMENT
Build the x_10513 plan (`plan10513.py` gives control seeds + detach set, already handles the
`X = T + p*t` choice that makes the pin multiplier divide) and solve its endgame EXACTLY instead of
greedily: enumerate the 3 atoms whose equations lie inside the 7-equation cluster, write the
7 x 3 integer system, and find the minimum-weight coset leader over the p-cosets of the detached
handle values.  My greedy closure4 only reaches 38,988 there because it cannot repair a688,
a19299 and the 1-equation shadow atoms (a16509, a39553, a41520, a41532) that contain the detached
variables — those need the exact solve, not greedy repair.

## Commands
python3 checker.py agentC_work/BEST_39013.json      -> 39013/39033
python3 agentC_work/DECISIVE.py                     -> the three-condition test
python3 agentC_work/carry2.py                       -> the P1 = P2 impossibility
python3 agentC_work/mincost.py ; cluster.py         -> the cost model
python3 agentC_work/run10513.py 10513 4             -> the best-cluster construction
