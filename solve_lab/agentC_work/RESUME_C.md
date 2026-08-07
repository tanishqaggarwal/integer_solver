# Agent C — RESUME (automated reasoning angle).  Updated after the decompilation.

## Verified baseline: best/new_instance_partial_39026.json = 39026/39033 (CONFIRMED by checker).
Installed (were absent): z3-solver, python-sat, cvc5, python-flint, ortools, sympy, numpy.
Rebuild caches first: `cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`

## THE INSTANCE, DECOMPILED (my main contribution; all in solve_lab/agentC_work/)
1. Break the 900 duplicate-equality 2-SCCs -> the circuit is a pure DAG with **8,173 free inputs**.
   Forward-evaluating from free inputs = 0: **0 broken gates, only 6 nonzero checks, score 39,005.**
2. Those 6 = three conditions: `x_18956 = K1 (mod p)`, `x_24468 = K2 (mod p)`, `OR(s1,s2)=1`.
   s1 = x_7715 = OR of 256 leaves, s2 = x_34554 = OR of 128 leaves.  In branch (1,1) the outputs are
   the FREE inputs x_22162, x_30213 -> set them to K2, K1.
3. `curve2.py/curve3.py/order.py`: the residual is exactly EC point addition,
   `y^2 = x^3 + a2 x^2 + a4 x + a6` (constants in agentC_work/curve.json), j = 0, short form A=0,
   **group order = the secp256k1 prime order n** ([n]G = O verified).  All 256 free leaf bits carry a
   pinned curve point and (`chain.py`) they are a **DOUBLING CHAIN P_i = 2^i G**.  The tree sums the
   selected leaf points; the root addition check binds `P1 + P2 = (x_22162,x_30213)`.
4. Two exits: (a) `P1+P2 = Q` = ECDLP on secp256k1 — dead; instance literals tested as the dlog
   (`consts.py`, 2,800 candidates) — no hit.  (b) **P1 = P2** makes A,B vanish identically and frees
   the root output.  P1=P2 needs kA = kB (mod n) with disjoint bit supports (178 on s1, 78 on s2),
   forcing kA-kB = +-n; both carry chains overflow (`carry2.py`) => **exactly impossible for free.**
5. `best_analyze.py`: the 39,026 deliverable FAKES P1 = P2 by overriding leaf-point pins; its residue
   sits on 7 handle-definition atoms {22229,22230,35758..35762} in 12 equations, 5 of which cancel.
6. `pairsweep.py` reproduces that construction for ANY of the 178x78 = 13,884 bit pairs
   (force P_u = P_w by overriding one leaf's pinned coordinate free-vars).  Pair (24601,2081) with
   closure2 gives 39,001; the gap to 39,026 is that closure2 cannot repair through a p-handle.
   `close3.py` adds that (solve for a gate output, then REALIZE it down the definer DAG).

## NEXT EXPERIMENT
Sweep `pairsweep.py sweep <shard> <n>` with closure3 over the 13,884 bit pairs (~2 s each).
Prior sessions only ever used one pair; the bit pair changes which handle atoms carry the residue and
therefore the equation count.  Also: choose the coordinate representative `C_w + p*t` with t making the
pin's small multiplier divide, which zeroes the pin atom and moves the residue purely onto the handle.

## Artifacts:  agentC_work/c3_*.json, PS_*.json, close_out.json (39,013), LOG.md
Best of mine so far: 39,013.  Anything >= 39,026 -> agentC_work/BEST_<score>.json immediately.
