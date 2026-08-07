# Agent C log — automated reasoning angle

## Step 0. Baseline verified
`python3 checker.py best/new_instance_partial_39026.json` -> satisfied 39026/39033, failing
[12231, 12270, 12350, 14584, 18673, 22044, 29125].  CONFIRMED.

Environment had NO solvers at all.  Installed: z3-solver 5.0.0, python-sat 1.9, cvc5 1.3.4,
python-flint 0.9.0, ortools 9.15, sympy 1.14, numpy.  Network works (pip).

## Step 1. Rebuilt s9 caches (atomize/poly/gates/fwd) — 57 s.  42,267 atoms, 39,033 eqs.

## Step 2. INDEPENDENT re-derivation of the circuit (agentC_work/supp2.py, supp3.py, fwd.py)
* greedy gate orientation leaves 1,800 vars "cyclic"; ALL 900 nontrivial SCCs are size 2 and are
  duplicated equalities `x_a = x_b`.  Breaking them gives a pure DAG with **8,173 free inputs**
  (not 7,273) and 30,575 gates, all with output coefficient exactly 1.
* **Forward evaluation from free inputs = 0 gives 0 division failures, 0 nonzero gate atoms,
  and only SIX nonzero check atoms.**  score 39,005 / 28 failing equations.

## Step 3. The six reduce to THREE independent conditions
```
a688   = 8863713*(x_18956 - K1) - x_14257      K1 = 1257873147476011081160397251633617631165504656759811518388115168273279192288235977446356 26
a40608 = a688^2                                 (redundant)
a1618  = x_24468 - K2 - x_32989                 K2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
a23000 = x_9274 - (x_29237 - x_23134)  = 1      x_9274 = x_2300 = 1 (pinned), x_29237-x_23134 = OR(x_7715,x_34554)
a39067, a41211                                  (both are a23000 times a constant; redundant)
```
* `x_23917 = x_22399 = p` exactly (p = 2^256-2^32-977).  `x_14257 = p*x_7497` with **x_7497 free
  and appearing in exactly one atom** => a688 = 0  <=>  **x_18956 = K1 (mod p)**.
* `x_32989 = p*x_11436` similarly => a1618 = 0 <=> **x_24468 = K2 (mod p)**.
* `x_18956 = x_37892 + p*h`, and `x_37892` is a 3-way MUX:
  `x_37892 = s1*(1-s2)*x_16742 + s2*(1-s1)*x_24908 + s1*s2*x_30213`, s1 = x_7715, s2 = x_34554.
* `x_7715 = OR of 256 leaves`, `x_34554 = OR of 128 leaves` (ortree2.py).  Condition C = at
  least one of those 384 leaves is 1.  Many leaves are FREE inputs; the rest are pinned to 0.
* `x_24468 = x_13913 + x_38045 + 12354891*p*h'`  =>  condition B is `x_13913+x_38045 = K2 (mod p)`.

## Step 4. Greedy topological closure (agentC_work/close.py) — 2.3 s to 39,013
Seeds {x_542:1, x_91:1, x_22162:K2, x_30213:K1} then repeatedly: for every nonzero check, pick the
latest-in-topological-order FREE variable occurring linearly and solve for it.  Rounds:
38,999 -> 38,978 -> 38,988 -> 38,986 -> 39,000 -> **39,013**, then STALL with 5 nonzero checks:
a19297 `x_15298*x_11150 + x_4007`, a19299 `x_15298*x_25739 - 6672769*x_29804`,
a30984 `537773*(x_15298*x_37758) - x_35605`, a36185, a40812.
Since x_15298 = s1*s2 = 1 this demands `x_11150 = x_25739 = x_37758 = 0 (mod p)` — i.e. exactly the
residual prior session 12 reported at §131.  Independently reproduced in 2.3 s from a cold start.
