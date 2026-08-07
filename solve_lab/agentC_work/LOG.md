# Agent C log — integer-polynomial analysis of EQUATIONS.txt

All statements below are about the raw integer and polynomial content of the file: atoms
(polynomials that must vanish), equations (integer linear combinations of atoms, some squared),
explicit integer literals, and congruences modulo the literal `p = 2^256 - 2^32 - 977` which
occurs in the file as the value of 220 variables.

## Step 0. Baseline verified
`python3 checker.py best/new_instance_partial_39026.json` -> satisfied 39026/39033, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.  CONFIRMED independently.
Environment had no solvers; installed z3-solver, python-sat, cvc5, python-flint, ortools, sympy,
numpy.  Caches rebuilt with `cd s9 && python3 atomize.py && poly.py && gates.py && fwd.py` (57 s):
**42,267 atoms, 39,033 equations, 38,748 variables.**

## Step 1. The instance is a triangular polynomial system (agentC_work/supp2.py, supp3.py, fwd.py)
* The greedy gate orientation leaves 1,800 variables "cyclic"; **all 900 nontrivial SCCs have size 2
  and are duplicated equalities `x_a = x_b`**.  Breaking them gives a strict partial order:
  **8,173 free variables** and 30,575 variables each defined by one atom, every one of those atoms
  linear in its defining variable with coefficient exactly 1 (so forward substitution never divides).
* Substituting free variables = 0 and evaluating forward: **0 of the 30,575 defining atoms is
  violated, and only SIX of the 10,792 remaining atoms are nonzero.**  Score 39,005.

## Step 2. The six reduce to three independent conditions
```
a688   = 8863713*(x_18956 - K1) - x_14257          a40608 = a688^2                (redundant)
a1618  = x_24468 - K2 - x_32989
a23000 = x_9274 - (x_29237 - x_23134) = 1          a39067, a41211                 (multiples of a23000)
K1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
```
* `x_23917 = x_22399 = p` exactly, and `x_14257 = p*x_7497` with `x_7497` free and occurring in
  exactly one atom.  Hence **a688 = 0 <=> x_18956 = K1 (mod p)**.  Likewise
  `x_32989 = p*x_11436`, so **a1618 = 0 <=> x_24468 = K2 (mod p)**.
* `x_9274 = x_2300 = 1` is pinned by a literal, and `x_29237 - x_23134 = s1 + s2 - s1*s2` where
  `s1 = x_7715`, `s2 = x_34554`.  Each of `s1`, `s2` is an iterated `u + w - u*w` over a binary tree
  of 256 and 128 leaves respectively (agentC_work/ortree2.py); **256 of those 384 leaves are free
  variables, the other 128 are pinned to 0**.  So a23000 = 0 asks that at least one of the 256 free
  leaf variables be 1.
* `x_18956` and `x_24468` are three-way selector combinations:
  `x_18956 = s1(1-s2)*x_16742 + s2(1-s1)*x_24908 + s1*s2*x_30213 + p*(free)` and the same shape for
  `x_24468` over `x_12186, x_14853, x_22162`.  **`x_16742, x_30213, x_22162, x_14853` are free
  variables**, so in the branch `s1 = s2 = 1` both conditions are met by setting
  `x_22162 = K2, x_30213 = K1` — verified EXACTLY OVER Z, not merely mod p.

## Step 3. Closure engines and the verified partial
`close.py` (greedy topological repair), `close2.py` (repairs re-solvable each round), `close3.py`
(solve an atom for a defined variable, then realise that value down the definition DAG),
`close4.py` (adds frame DETACHMENT: drop a variable's defining atom so the variable becomes a
control input and its old definer becomes a violated check).  From the seeds above, closure reaches
**39,013 in 2.3 s from cold** — `agentC_work/BEST_39013.json`, **checker-verified**
(20 failing).  Residual: `a19297 = x_15298*x_11150 + x_4007`,
`a19299 = x_15298*x_25739 - 6672769*x_29804`, `a30984 = 537773*(x_15298*x_37758) - x_35605`,
plus `a36185`, `a40812` (one equation each).  Since `x_15298 = s1*s2 = 1` these demand
`x_11150 = x_25739 = x_37758 = 0 (mod p)`.

## Step 4. The residual, in closed polynomial form
`x_11150, x_25739, x_37758` are rank-2 combinations of two variables `x_35389` and `x_6671`, whose
values I measured against candidate polynomials over random probes.  Writing
`u1 = x_12186, v1 = x_16742, u2 = x_14853, v2 = x_24908, u3 = x_22162, v3 = x_30213`:
```
x_6671  = (v3 + v1)*(u2 - u1) - (u1 - u3)*(v2 - v1)                 exact on every probe
x_35389 = (u2 - u1)^2 * (u3 + u1 + u2 + c) - (v2 - v1)^2            exact on every probe
c = 97553848499418123410591666447050222001188385549510401465815187079080512838891
```
`u1` and `v2` are constant under every probe of the free variables `v1`, `u2`.  The identities are
recorded as measured polynomial facts about the file; no interpretation is attached to them.
`x_35389 = x_6671 = 0` holds identically when `u1 = u2` and `v1 = v2`, which makes the residual
vanish with `u3, v3` unconstrained — that is what the 39,026 deliverable exploits.

## Step 5. Why `u1 = u2, v1 = v2` is not free (agentC_work/site.py, bitcost.py, truecost.py)
* `v1 = x_16742` is pinned by `a26731 = 6788513*(x_16742 - x_19083) - x_9254` with
  `x_9254 = p*(free)`, so `x_16742 = x_19083 (mod p)`; likewise
  `a29539 = 12846437*(x_14853 - x_1308) - x_29967` pins `u2` to `x_1308 (mod p)`.
* Each of the 256 free leaf variables `b` carries conditional pins `b*(X - C) - m*H` with `C` an
  explicit literal, `m` a small literal and **`H = p*(free)` in 512 of 512 cases**.  So setting
  `b = 1` pins `X` modulo `p` to a literal of the file.  Moving those pinned values requires
  detaching a handle-definition atom, which violates it.
* `best_analyze.py`: the 39,026 deliverable does exactly that.  Its violated atoms are
  `{22229, 22230, 35758, 35759, 35760, 35761, 35762}`, living in 12 equations of which 5 vanish by
  cancellation, giving its 7 failures.

## Step 6. The x_10513 cluster: predicted 39,027, REFUTED by exact computation
`bitcost.py`/`truecost.py` rank the 256 leaf variables by the union of equations their two handle
atoms touch; the minimum is **7 at x_10513** (a8427: 7 equations, a8429: 5, union 7) against 11 for
the deliverable's x_24601.  With the balance law `failing = |E| - n + c` and `n = 8` atoms whose
whole equation set lies inside `E`, that predicted 6 failures (39,027).

`exact10513.py` builds the 12 x 8 incidence matrix explicitly:
```
eq748   [ 25  13  14 | 0 0 0 0 0]      eq1666  sq [0 0 0 | 1 0 0 0 0]
eq1785  [ -3 -13 -14 | 0 0 0 0 0]      eq12466 sq [0 0 0 | 0 1 0 0 0]
eq2629  [ 38  27 -22 | 0 0 0 0 0]      eq26941 sq [0 0 0 | 0 0 1 0 0]
eq3676  [-27 -36   0 | 0 0 0 0 0]      eq30004 sq [0 0 0 | 0 0 0 1 0]
eq5692  [  1 -29 -28 | 0 0 0 0 0]      eq30122 sq [0 0 0 | 0 0 0 0 1]
eq5717  [ 20 -37 -18 | 0 0 0 0 0]      eq20538    [30 0 0 | 0 0 0 0 0]
```
**rank = 8 = number of atoms, kernel ZERO.**  The five atoms counted as compensators
(a16509, a39553, a41277, a41520, a41532) are **SHADOWS** — each a fixed linear combination of
a8427/a8428/a8429 sitting alone in its own squared equation.  They inflate `n` without adding one
degree of freedom, so the balance law does not apply here.  True minimum for the cluster is
**>= 11 failures**; the measured construction gives 38,988.  **The 39,027 prediction is refuted.**
Note **eq20538 = 30*a8427 alone** — a single-atom equation forcing that atom to zero, structurally
identical to **eq29125 = a22230** in the deliverable's cluster.  The file places one such guard in
every cheap cluster found so far.

## Step 7. Classifier-free placement bound (agentC_work/minweight.py) — the live scan
My earlier `globalscan.py` priced clusters by counting "settable" atoms of the shape
`out = p*(free)`, which scores the deliverable's own cluster at 9 when it achieves 7, so its
no-better-cluster verdict was suggestive only.  The fix removes the classifier entirely:
for a cluster (E, S) the number of failures is `|{e in E : M_e . v != 0}|`, and minimising that over
ALL nonzero rational `v` is a minimum-weight-codeword problem in the code generated by the columns
of `M`.  It ignores which atoms are settable and ignores integrality, so it is a RELAXATION and the
true cost is `>=` it.  Any cluster whose min weight is `>= 7` therefore cannot beat 39,026, with no
settability judgement involved.  Minimum weight is attained on a `v` annihilating `|S|-1`
independent rows, so those are enumerated exactly.  Calibration target: the deliverable's own
cluster must return 7.

## Step 8. CALIBRATION RESULT — the requested classifier fix is IMPOSSIBLE as posed
```
DELIVERABLE cluster: |E|=12 |S|=8  ->  classifier-free min weight = 5   (observed failing = 7)
x_10513     cluster: |E|=7  |S|=3  ->  classifier-free min weight = 4   (true cost >= 11)
```
The task was to fix the "settable atom" classifier so it reproduces the deliverable's cost of 7 and
then rescan.  The classifier-free relaxation is the strongest purely structural pricing available
(it drops settability AND integrality, so it is a sound lower bound), and it returns **5, not 7**.

**The gap of 2 is exactly the congruence term, and it is not structural.**  Every handle in the file
has the form `p*(free variable)` — verified for 512 of 512 leaf pins and for the two coordinate pins
a26731/a29539 — so each violated atom's value is confined to a FIXED residue class mod p determined
by which literal constants the construction lands on.  A rational kernel direction of `M` is
therefore realisable only if it is compatible with those residues, which costs the 2 equations.
Since the residue vector depends on the construction and not on the equation-atom incidence matrix,
**no function of the incidence structure alone can return 7.**

Consequences, stated plainly:
* My earlier `globalscan.py` under-counted; the corrected pricing over-counts freedom instead.  The
  no-better-cluster verdict is therefore **still suggestive, not established** — I could not upgrade
  it, and I am recording that rather than presenting the relaxation as if it settled the question.
* The relaxation is also very loose (4 against a measured >= 11 on x_10513), so it is not useful for
  ranking either.  Structural pricing of defect placement is a dead instrument.
* The version that would work is residue-aware: compute `v0 = atom values mod p` from a real
  construction, keep only subsets `T` with `M_T v0 = 0 (mod p)`, then solve over the p-lattice.
  That equals 7 on the deliverable by construction.  It costs one construction per cluster, so it
  must be aimed at a few dozen candidates, not scanned over 3,349.
