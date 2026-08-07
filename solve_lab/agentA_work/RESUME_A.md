# Agent A — RESUME (exact integer linear algebra / lattice angle)

## Best verified score: 39,026 (baseline, re-verified myself)
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json` -> 39026/39033,
failing [12231,12270,12350,14584,18673,22044,29125].  No improvement of my own yet.

## Setup needed to re-enter (caches were MISSING from the repo)
```
cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py
cd ../agentA_work && python3 probe5.py     # dumps the 33 region atoms
```

## What I established (re-derived from the file, not inherited)
1. The residual region is EXACTLY 33 atoms in 39 equations (12 = E, 27 outside, all
   currently satisfied). `region.py` computes the closure levels:
   L0 33 atoms/11 knobs/39 eqs; L1 471/341/1105; L2 5807/4045/14225.
2. The 33 atoms are five parallel CHAINS of the shape (p-handle atom, check atom):
   5766-5772 | 10935-10938 | 19087-19092 | 22229-22235 | 35754-35762  (see probe5.py output).
3. **NEW knobs prior work did not use**: `region.py` says 11 variables have ZERO atoms
   outside the 33 (prior generator list had 9, and only over the seven):
   642, 1329, **1613**, **1844**, 9413, 10903, 17325, **21574**, **29305**, 29854, 31864.
   A further ~13 vars have exactly ONE outside atom, and for eight of them that atom is
   a37887 which lives in exactly ONE equation (8680): 950, 6947, 9629, 15120, 23754,
   28730, 33168, 35619.  x2892/x6090/x28355 -> a41906 only.
4. `ahandles.py`: 1,562 atoms carry a private variable (occurs in no other atom);
   **326 have granularity 1** (fully free atom values).  s10/handles.py counted only
   FREE INPUTS and so reported 1,249/all-p-quantised — my census is strictly larger.
   None of the 326 is in E, so this does not by itself break the residual.
5. Exact restatement of the obstruction: reachable alpha (7 atom values) =
   {alpha1 = 0 mod p fixed residue K2} and {alpha0 + 7376877*alpha6 = C0 mod p}.
   12 rows, rank 7 -> all 12 need alpha = 0 mod p -> impossible while C0,K2 != 0.
   Breaking ONE congruence buys exactly 1 equation (39,027); breaking BOTH gives 39,033.

## Next experiment (in progress)
Exact max-satisfy integer program on the L0 region ENLARGED by a37887/a41906:
40 equations x 33 atom values, knob lattice from the 19-24 knobs above; enumerate
dropped-subsets of size <= 6 and test integer solvability by HNF.  Script: `maxsat.py`.
