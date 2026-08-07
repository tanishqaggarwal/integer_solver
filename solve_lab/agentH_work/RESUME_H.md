# Agent H — RESUME (decomposition angle).  STOPPED by coordinator.

## Best verified score
39,026 — the inherited deliverable `solve_lab/best/new_instance_partial_39026.json`, re-verified by me
(`satisfied 39026/39033`, failing [12231,12270,12350,14584,18673,22044,29125]).  I did NOT beat it.
My own best built from scratch: **39,018** (`scan1.py`, any single bit; 4 nonzero atoms
{30980,30982,36185,40812}); also 39,014 saved at `E_39014_542_1438.json`.  Nothing of mine >= 39,026.

## THE FRAME (main durable result; see DECOMPOSITION.json)
Orient every atom of syntactic form `x_t - rest` as the definition `x_t := rest`.
=> gate DAG is **ACYCLIC**, covers all 38,748 vars, **8,747 free inputs**, 30,001 defined,
12,266 check atoms.  (Prior sessions' frame: 7,273 free inputs, 1,800 vars in cycles.)
Forward eval from ALL-ZERO free inputs scores **39,005** with only **5** nonzero check atoms.
Files: model.py -> fwd2.py -> support.py -> ev.py (exact) / fast.py (incremental, 1.4 ms/move).

## DECOMPOSITION MEASUREMENTS (DECOMPOSITION.json)
- eq-var bipartite graph: 1 component.  atom-eq graph: 1 giant + 3,234 singletons.
- Free-input hypergraph (hyperedge = equation): **1 component**, 8,747 vars, 18,248 equations.
  20,785 equations are identically satisfied by forward eval.
- Var free-support mean 2.5, max 259; equation free-support mean 5.9, max 282.
- Residual closure = 6,007 free inputs / 6,026 checks / 9,244 equations.  **No small separator.**

## WHAT THE INSTANCE IS (decompiled, mine, independent)
Residual = exactly 3 conditions: `OR(256 bits) = 1`, `x_37892 = C1 (mod p)`, `x_13682 = C2 (mod p)`.
The circuit is a binary MUX TREE over 256 bits.  Each bit b carries a leaf point P_b = its two
load-pin constants.  At a node where both children fire the free drivers (X3,Y3) are UNPINNED
(gate `1-L*R`) but three checks of rank 2 force the chord identities
`(X1+X2+X3+K)(X2-X1)^2 = (Y2-Y1)^2`, `(Y3+Y1)(X2-X1) = (Y2-Y1)(X1-X3)` mod p.
=> the tree computes an **elliptic-curve multi-addition**.
Curve (shifted x = X + K/3): y^2 = x^3 + B, B = 6401953368003087640844319876221082905875170063455428218598732582039359852479 4,
p = 2^256-2^32-977, group order = secp256k1's n, PRIME.  All 256 leaf points and the target are on it.
**255 of 256 points satisfy 2*P_i = P_j: one doubling chain of length 256 starting at bit x_2779.**
So P_i = 2^i G and the instance is satisfiable iff `sum_{i in S} 2^i * G = P*`, i.e. **S is the binary
expansion of the discrete logarithm of P* base G on a 256-bit prime-order curve.**
No subset of size <= 4 hits P* (exhaustive).

## Re-entry
cd solve_lab/agentH_work
python3 model.py; python3 fwd2.py; python3 support.py   # rebuild caches (~1 min)
python3 close2.py <u-bit> <w-bit>   # construct+close a branch; python3 scan1.py  # 39,018 states
python3 ec.py                       # curve/points; expo.json = bit -> exponent i in P_i = 2^i G

## Single next experiment
Test whether the discrete log k is weak before conceding: (a) BSGS for k < 2^44,
(b) meet-in-the-middle for Hamming weight <= 6 over the 256 exponents, (c) runs of consecutive
set bits (k = 2^a(2^m-1)), (d) k congruent to a small multiple of a subgroup-free structure.
If all fail, 39,033 requires a 256-bit ECDLP and the correct deliverable stays 39,026.
UNVERIFIED: my 2-bit test of "delivered point == P_b1 + P_b2" did not confirm, because the closer
left the three EC checks nonzero (it had the drivers frozen).  Re-run with X3,Y3 unfrozen.
