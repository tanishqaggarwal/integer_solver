# Agent H — RESUME (decomposition angle)

## Best verified score
39,026 (inherited `solve_lab/best/new_instance_partial_39026.json`, re-verified by me).
My own best from scratch: 39,014 (`E_39014_542_1438.json`) — not yet competitive, but built from a
*complete structural decompilation* of the instance (below), reached in one pass with no search.

## What I established INDEPENDENTLY (all from my own parse)
1. `model.py` -> 42,267 atoms / 39,033 eqs.  `fwd2.py`: orienting every atom of syntactic form
   `x_t - rest` as a DEFINITION gives an ACYCLIC gate DAG covering all 38,748 vars with
   **8,747 free inputs** and 12,266 check atoms.  (Prior sessions' frame had 7,273 free inputs and
   1,800 vars in cycles — my frame is strictly cleaner and has no cycles.)
2. Forward eval from ALL-ZERO free inputs scores **39,005** with only **5 nonzero check atoms**
   (`ev.py`, verified by solve_lab/checker.py on `allzero_fwd.json`).
3. Variable free-input supports are tiny (mean 2.5, max 259); 20,785 of the 39,033 equations are
   identically satisfied by forward eval.  Reduced problem = 18,248 eqs over 8,747 free inputs,
   ONE connected component (no free graph decomposition; `decomp.py`).
4. **The entire residual is 3 conditions**: `x_9274 = OR(all 256 bits) = 1`,
   `x_37892 ≡ C1 (mod p)`, `x_13682 ≡ C2 (mod p)`, p = 2^256-2^32-977,
   C1 = 12578731474760110811603972516336176311655046567598115183881151682732791922882359774463562**6**,
   C2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002.
5. **DECOMPILATION**: the instance is a binary MUX TREE over 256 boolean bits computing an
   ELLIPTIC-CURVE MULTI-ADDITION.  Each bit b carries a leaf point P_b = (H1_b, H2_b) (its two
   load-pin constants).  At every tree node with both children firing, the node's free drivers
   (X3,Y3) are UNPINNED (gate is `1 - L*R`) but three checks force the chord identities
   `(X1+X2+X3+K)(X2-X1)^2 = (Y2-Y1)^2` and `(Y3+Y1)(X2-X1) = (Y2-Y1)(X1-X3)` mod p,
   K = 97553848499418123410591666447050222001188385549510401465815187079080512838891.
   Hence the delivered value = **EC sum of the points of all set bits**, and the instance is
   satisfiable iff **SUM_{b in S} P_b = P* = (C2, C1)** in shifted coords x = X + K/3.
6. `close2.py`/`close3.py`: a bottom-up cascade closer (smallest-atom-support first, freeze each
   assigned var) closes the whole pin tree in ~32 exact assignments, no search.

## Re-entry
cd /home/user/integer_solver/solve_lab/agentH_work
python3 model.py; python3 fwd2.py; python3 support.py; python3 ev.py   # rebuild caches
python3 close2.py <u-bit> <w-bit>     # construct + close a 2-bit branch

## Next experiment (highest value)
Extract the 256 leaf points, verify they lie on one curve, then solve the EC subset-sum
SUM_{b in S} P_b = P*: single bits, pairs, triples, then meet-in-the-middle.
