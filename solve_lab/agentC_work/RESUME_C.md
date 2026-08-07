# Agent C — RESUME (automated reasoning angle)

## Verified baseline
`python3 checker.py best/new_instance_partial_39026.json` -> 39026/39033. CONFIRMED.
Installed (were absent): z3-solver, python-sat, cvc5, python-flint, ortools, sympy, numpy.
Rebuild caches first: `cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`

## My independent re-derivation (much simpler than prior sessions' narrative)
Scripts in `solve_lab/agentC_work/`: supp2.py (SCC), supp3.py (support/DAG), lib2.py (fast fwd eval,
24 ms), ort.py (OR-tree leaf extractor), cone.py/trace.py/var.py (inspectors), try1.py (constructions).
* All 900 nontrivial SCCs are size-2 duplicate equalities -> circuit IS a DAG, **8,173 free inputs**.
* Forward-evaluating from **free inputs all zero** gives 0 broken gates and only **6 nonzero checks**
  (score 39,005).  The 6 reduce to 3 independent conditions:
  - A: `x_18956 = K1 (mod p)`  (a688; x_14257 = p*x_7497, x_7497 free & solo)
  - B: `x_24468 = K2 (mod p)`  (a1618; x_32989 = p*x_11436)
  - C: `OR(x_7715, x_34554) = 1` (a23000). s1=x_7715=OR of 256 free/pinned leaves, s2=x_34554=OR of 128.
* `x_18956 = MUX(s1,s2)[x_16742, x_24908, x_30213] (mod p)`, `x_24468 = MUX(s1,s2)[x_12186,x_14853,x_22162]`.
* **x_16742, x_30213, x_22162, x_14853 are FREE inputs.**  Branch (s1,s2)=(1,1) selects (x_22162,x_30213),
  both free -> set x_22162=K2, x_30213=K1.
* CONSTRUCTED (try1.py): seeds {x_542:1, x_91:1, x_22162:K2, x_30213:K1} -> **A and B both TRUE**,
  score 38,999, only 6 nonzero checks left, and all six are the *conditional pins* of the two bits:
  `b*(X-C) - m*h`.  The X's (x_13153, x_20386, x_35344, x_23210) are FREE inputs -> set them to C.

## NEXT EXPERIMENT (running now)
Set those 4 free vars to their pin constants on top of the above.  If the handles land at 0 this is a
FULL SOLVE. If not, repeat the "activate bit -> satisfy its 2 pins" closure until fixpoint.

## Artifacts
Best so far is still `best/new_instance_partial_39026.json` (not mine). Anything >=39,026 of mine goes to
`agentC_work/BEST_<score>.json` immediately.
