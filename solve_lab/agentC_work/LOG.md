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

## Step 5. THE INSTANCE DECOMPILED (agentC_work/curve2.py, curve3.py, order.py)
* The residual after closure is the point-addition pair
  `A = (x2-x1)^2*(x3+x1+x2+a2) - (y2-y1)^2`,  `B = (y3+y1)(x2-x1) - (x1-x3)(y2-y1)`,
  measured to match `x_35389`/`x_6671` digit for digit, with
  `x1=x_12186, y1=x_16742, x2=x_14853, y2=x_24908, x3=x_22162, y3=x_30213` and
  `a2 = 97553848499418123410591666447050222001188385549510401465815187079080512838891`.
* Fitting `y^2 = x^3 + a2 x^2 + a4 x + a6` from P1,P2 and testing Q: **Q lies on it** (nontrivial),
  a4 = 114170008767671698752186727197936107864370654164657728518655355473804451402762
  a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490
* Short form X = x + a2/3 gives **A_short = 0**, B_short = 6401953368003087640844319876221082905875170063455428218598732582039359852479 4
  -> j = 0, and `[n_secp]G = O` with n_secp = the secp256k1 group order (prime).
  **The curve is isomorphic to secp256k1; the group is the prime-order secp256k1 group.**
* All **256 free leaf bits** carry a pinned point on that curve (curve3.py: 256/256), and
  `P1`,`P2` in the solved state ARE leaf points.  The selector tree is a binary accumulation tree:
  root value = SUM of the selected leaf points; the A/B checks bind only at nodes where both
  children are active.
* **Therefore the instance = "find a subset of 256 given secp256k1-isomorphic points summing to Q".**
  s1-side has 178 free bits, s2-side 78.

## Step 7. The cost model, computed exactly (agentC_work/mincost.py, cluster.py, close4.py)
* Every free leaf bit b carries two conditional pins `b*(X - C) - m*H` whose handle H is
  `p * (free var)`.  So when b = 1, X is pinned mod p and the leaf point is fixed mod p.
  Faking `P1 = P2` therefore requires DETACHING the handle-definition atom(s) (frame change),
  which is what `close4.py` implements (`detach` set: the definer is dropped, the variable becomes
  a controllable input, and the definer atom becomes a broken check).
* `mincost.py`: over all 256 bits, `min |E(hx) u E(hy)| = 7` at **bit x_10513**
  (handles a8427 = 7 eqs, a8429 = 5 eqs; union 7).  The 39,026 deliverable's bit x_24601 costs 11.
  Cross plans (override u's x and w's y) all cost >= 12, so the same-bit plan is optimal.
* `cluster.py` (balance law `failing = |E| - n + c`, n = atoms whose whole equation set is inside E):
  - deliverable cluster {22229,22230,35758..35762}: |E| = 12, n = 8, |E| - n = **4**
  - **bit x_10513 cluster {a8427,a8429}: |E| = 7, n = 3, |E| - n = 4**  <- same slack, half the size
  So the x_10513 cluster admits `failing = 7 - 3 + c`; with the deliverable's c = 2 this is **6**
  (score 39,027) and with c = 0 it is 4 (score 39,029).  That is the concrete route past 39,026.
* First construction on it scored only 38,989 because the greedy closure could not repair a688,
  a19299 and the 1-equation "shadow" atoms (a16509, a39553, ...) that contain the detached
  variables.  Patched close4 to forbid touching the 220 p-wires during realize().

## Step 8. DECISIVE CROSS-CHECK (agentC_work/DECISIVE.py) — answers to the coordinator
Q1: does branch (1,1) with x_22162 = K2, x_30213 = K1 close the system?
    NO — but it satisfies all three top conditions EXACTLY OVER Z (a688 = a1618 = a23000 = 0).
    The only checks left nonzero are the two activated bits' own conditional pins; closing those
    (close2.py) gives **39,013, checker-verified** (agentC_work/BEST_39013.json).
Q2: is the remaining breakage exactly the point-addition law?  YES, exactly:
    a19297 `x_15298*x_11150 + x_4007`, a19299 `x_15298*x_25739 - 6672769*x_29804`,
    a30984 `537773*(x_15298*x_37758) - x_35605`, plus a36185 and a40812 (1 eq each);
    x_11150/x_25739/x_37758 are rank-2 in A = x_35389 and B = x_6671, and I measured
    B = (y3+y1)(x2-x1) - (x1-x3)(y2-y1) matching digit for digit over random probes, and
    A = (x2-x1)^2 (x3+x1+x2+a2) - (y2-y1)^2 with a2 constant across probes.
Refinement vs agent I: the group order is EXACTLY n_secp ([n]G = O verified), so the curve is
ISOMORPHIC to secp256k1 itself, not a different-order sextic twist.  j = 0, short form A = 0.
Sharpening: the DLP is not the only door.  `P1 = P2` makes A and B vanish IDENTICALLY and frees the
root output, needing no dlog.  I closed that door exactly (carry2.py): P1 = P2 requires
kA - kB = +-n with disjoint bit supports and BOTH deterministic carry chains overflow.  So the
ECDLP is forced only after this second, purely combinatorial door is proved shut.

## Step 9. PRIORITY-1 RESULT: the 39,027 prediction is REFUTED BY MY OWN EXACT COMPUTATION
`exact10513.py` builds the x_10513 cluster explicitly: |E| = 12 equations, |S| = 8 atoms inside.
The 12 x 8 coefficient matrix is:

    eq748   [ 25  13  14  0 0 0 0 0]      eq1666  sq [0 0 0 1 0 0 0 0]
    eq1785  [ -3 -13 -14  0 0 0 0 0]      eq12466 sq [0 0 0 0 1 0 0 0]
    eq2629  [ 38  27 -22  0 0 0 0 0]      eq26941 sq [0 0 0 0 0 1 0 0]
    eq3676  [-27 -36   0  0 0 0 0 0]      eq30004 sq [0 0 0 0 0 0 1 0]
    eq5692  [  1 -29 -28  0 0 0 0 0]      eq30122 sq [0 0 0 0 0 0 0 1]
    eq5717  [ 20 -37 -18  0 0 0 0 0]      eq20538 [30 0 0 0 0 0 0 0]

**rank(M) = 8 = |S|, so the kernel is ZERO.**  The five atoms I had counted as free compensators
(a16509, a39553, a41277, a41520, a41532) are **SHADOWS** — each is a fixed linear combination of
a8427/a8428/a8429 sitting alone in its own square equation.  They inflate n without adding any
freedom.  The balance law `failing = |E| - n + c` therefore does not apply to this cluster:
the correct count is  failing = (equations whose linear form in the two detached values is
not identically zero) = 6 bilinear rows + eq20538 + 4 shadow squares = **>= 11**, and the mod-p
residues of a8427/a8429 are pinned by the pin construction so no row can be cancelled
(each cancellation needs a 1-in-p coincidence).  Measured construction agrees: 38,988.
Also note **eq20538 = 30*a8427 alone** — a single-atom equation forcing a8427 = 0, exactly the
same structural guard as eq29125 = a22230 in the deliverable's cluster.

`globalscan.py` over all 3,349 settable handle-definition atoms: the best single-seed cluster is
a8429 with |E| = 5 and rank 1 (|E| - rank = 4), but breaking a8429 alone cannot fake P1 = P2
(both coordinates must move).  My "settable" classifier scores the deliverable's own cluster at
|E| - rank = 9 while it actually achieves 7, so the classifier under-counts the deliverable's
freedom (its 7 atoms are settable through free inputs x_9118, x_8731, x_642, x_17325, x_1329,
x_10903, x_9413, not through the wire*free pattern).  **No cluster beats 12/7; 39,026 stands.**

## Step 10. CURVE RECONCILIATION (settles agent G vs agent I vs me)
From my curve.json alone, all three descriptions are the same curve:
* The circuit's law is `x3 = lambda^2 - a2 - x1 - x2`, i.e. GENERAL Weierstrass with a nonzero x^2
  coefficient.  Agent G's "extra constant K" **is** a2; it is not an anomaly, it is the a2 term.
* Testing `y^2 = x^3 + 7` is therefore the wrong form.  The real form is
  a2 = 97553848499418123410591666447050222001188385549510401465815187079080512838891
  a4 = 114170008767671698752186727197936107864370654164657728518655355473804451402762
  a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490
  and all 256 leaf points, P1, P2 and **Q** lie on it (Q on y^2=x^3+7: False; Q on the real form: True).
* Substituting X = x + a2/3 gives A_short = **0** exactly, so j = 0, and
  B_short = 64019533680030876408443198762210829058751700634554282185987325820393598524794.
* `[n_secp]G = O` with n_secp = 115792089237316195423570985008687907852837564279074904382605163141518161494337,
  and **B_short/7 is a sixth power mod p** (verified).  So the curve is the TRIVIAL sextic twist:
  F_p-isomorphic to secp256k1 by x -> u^2 x, y -> u^3 y with u^6 = B_short/7.
  Agent I is right that b differs, but a *different b with the same group order* is an isomorphic
  copy, not a nontrivial twist (a nontrivial sextic twist would carry one of the other five CM orders).
