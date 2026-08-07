# Agent E log

## Pivot notice (mid-session)
Coordinator directed: DROP generator/authorship reverse-engineering (no PRNG hunting, no
emission-order forensics, no "run the recovered generator forward").  Work the equations as
given, algebraically/structurally.  **Withdrawn:** nothing was built on emission-order or
PRNG forensics; no such artifact exists.  `parse3.py` decomposes each polynomial into its
written sub-expression summands — this is a syntactic decomposition of the given polynomial
(a property of the mathematical object), and constraint propagation over it is item 1 of the
PROMPT's attack portfolio, so it is retained.  All conclusions below are computed from the
equations themselves.

## 0. Baseline verified
`python3 checker.py best/new_instance_partial_39026.json` -> 39026/39033,
fails [12231, 12270, 12350, 14584, 18673, 22044, 29125].  CONFIRMED.

## 1. Independent parse (parse3.py -> model3.pkl)
39,033 equations; 9,710 are perfect squares `S*S`, 29,323 linear.
Each equation = outer_scalar * (integer-linear combination of *atoms*).
40,727 distinct atoms over 38,748 variables.

## 2. Atom-level reproduction of the deliverable (aeval.py)
At `best/new_instance_partial_39026.json` exactly **8 atoms are nonzero**, and evaluating the
39,033 equations purely from those atom values reproduces the failing set exactly.
=> the atom model is faithful.

## 3. Gate DAG (dag.py -> dag.pkl)
35,004 atoms are of the form `x_out - RHS(x...)` with `x_out` not in RHS (a definition);
30,383 variables get defined (25,762 once, 4,621 twice); 8,365 variables are never defined
("free").  The definition relation is ACYCLIC.

## 4. **KEY RESULT — propagation from free=0 leaves only THREE violated atoms** (prop2/prop3)
Seed all 8,365 free variables to 0 and run unit propagation over atoms (linear + quadratic
single-unknown solve).  Result: all 38,748 variables determined, 40,724 / 40,727 atoms exactly
zero.  Violated atoms:
  a747   : x_14257 - x_23917*x_7497        (cone 11 vars, 2 free inputs: x_7497, x_18956)
  a20212 : x_13913 - (x_608 + x_22978)     (cone 2,209 vars, 261 free inputs)
  a24403 : x_9274 - (x_29237 - x_23134)    (cone 2,182 vars, 256 free inputs; residual = 1)
Checker on that assignment: **38,998 / 39,033** (35 failing).
This is a structurally much simpler residual than the 8-atom deliverable state.

## 5. Decoding the residual (trace.py / loc.py / cone dumps)
`a24403` = `1 = (a+b) - a*b` with `a=x_7715`, `b=x_34554` booleans -> **OR(a,b)=1 forced**
(`x_9274 = x_2300 = 1` is pinned by a literal).  `a=OR-tree over 88+90 free bits`,
`b=OR-tree over 78 free bits` (pure alias/sum/prod OR gadgets, verified by shape census).
`a20212` is a 2-way MUX `x_13913 = a(1-b)x_12186 + b(1-a)x_14853 (+ ab-term x_22162)`;
`a20215` is `x_24530 = b(1-a)x_24908`.

## 6. The (1,1) branch closes all four core atoms EXACTLY
Activate one free bit in each OR-tree, then set the free vars
`x_22162 = x_13682`, `x_30213 = x_18956 - x_32237` (4-iteration fixpoint).
Bits (x_4279, x_26005): bad set becomes exactly the two bits' pin atoms
`[6668, 12606, 34497, 34498]`, 21 failing equations.  Also works for (1,0):
`x_12186 = x_13682`, `x_16742 = x_18956 - x_32237` -> bad `[4872,4877,<2 pins>]`.

## 7. Pin gadgets, and why single-move greedy stalls at 4 atoms
Each activated bit b carries `b*(free - K) = m*handle` and `pin = p*handle`
(p = the 256-bit literal 115792089237316195423570985008687907853269984665640564039457584007908834671663 that appears throughout the file).  Fixing one pin relocates the residual one link down a chain
(`runs/chain_A.log`: 34498->10372->35326, 34497->13311->34090->35321,
12606->19725->7128, 6668->209).  #bad atoms is invariant = 4 under single moves.

## 8. The repair is an exact LINEAR DIOPHANTINE system  (jclose2.py / jsolve3.py)
`fast.py` gives incremental downstream-only re-evaluation (0.08 s/probe, exact).
Closure of the map free-vars -> atom residuals around bad=[6668,12606,34497,34498]:
converges in 4 s to **4,008 variables x 2,996 atoms**, 1,810 (var,atom) nonlinear pairs.
- round<=1: 61 vars / 86 linear rows -> rank(A)=rank([A|b])=59, **rationally feasible,
  integer-INfeasible** (kernel 2).
- round<=2: 342 vars / 335 linear rows -> rank 297 = rank([A|b]), rationally feasible; HNF running.
=> the obstruction at small support is INTEGRALITY, not rank.  Growing the support is the attack.

## 9. RESUMED — the bit scan and the subset test (the decisive experiments)

### 9.1 Machinery fixed
`sparse.py` rewritten with a column index + heap pivot order and guards
(`maxcore`, entry-bit-length); `bitfeas2.py` closes and solves per bit at growing support;
`scanbits.py` sweeps a whole tree with a per-bit alarm.

### 9.2 Per-bit pin feasibility  (`runs/scanB.log`, `runs/recheckB.log`, `runs/scanA5.log`)
For a single activated bit the pin system is an exact linear Diophantine system.
**b-tree: 78/78 scanned -> 50 FEASIBLE, 28 infeasible.**  Every feasible one was *exactly*
verified: applying the solution leaves **zero** atoms outside the selector core
(`exact=(16, [])` -- 16 failing equations, all from the still-unresolved selector atoms).
Infeasible bits fail on a p-divisibility row, e.g. `x_30448 -> row 23194: rhs % -p != 0`,
or on an infeasible core.  A re-check of all 28 with the improved solver reproduced the
verdict.  a-tree scan running the same way.

### 9.3 The subset law (`runs/subsetsize.log`, `runs/iterpair.log`)
Iterated closure+solve from k activated bits, residual atoms after convergence:

| k | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| residual atoms | **0, 0, 0** | 3, 3, 4 | 7, 7, 7 | 6, 7, 9 | 9, 10, 10 |

k=1 closes completely.  Every k>=2 leaves whole **AND-triples** and the iteration is a fixed
point on them (`runs/iterpair.log`: it1 makes no progress).

### 9.4 What an AND-triple is  (decoded exactly, pair (x_1530, x_1603))
Two bits in different halves of the OR-tree make their LCA's AND gate fire:
`x_24195 = x_33953 * x_1250 = 1`.  Then three atoms appear, in one shape:

    a726 : x_24195 * x_19097            = 0            x_19097 = 5002401*U  + 15322661*V
    a722 : x_24195 * x_6635             = p * x_34496  x_6635  = 15944455*U +  4826103*V
    a724 : 7952523*(x_24195 * x_36280)  = p * x_3193   x_36280 = 14913407*U + 11707765*V

with `U = x_29210 = x_25848 - x_17317`, `V = x_8736 = x_18682 - x_28841`, and
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663.
The 2x2 determinant of the two mod-p forms is 15944455*11707765 - 4826103*14913407
= **114700293930154**, and `gcd(114700293930154, p) = 1`, so the two congruences have rank 2
and force **U = V = 0 (mod p)**; a726 then forces `5002401*U + 15322661*V = 0` over Z.
This is the SAME shape as the top-level (1,1) residual `x_15298 * {x_11150, x_25739, x_37758}`.
At the repaired pair state the only free handles on (U,V) are 23 further OR-tree selector
bits -- moving them fires more AND nodes.

**Conclusion: the pin conditions are independent per single bit but NOT across subsets.
Any on-set of size >= 2 binds through the AND (LCA) nodes it activates, each demanding a
rank-2 mod-p vanishing of the same two forms.**

### Framing note (coordinator direction change)
All curve/point/group vocabulary has been removed.  `p` below is simply the 256-bit integer
literal 115792089237316195423570985008687907853269984665640564039457584007908834671663 that
occurs verbatim in EQUATIONS.txt; every statement here is about integer congruences and
integer-linear relations among the polynomials in the file, nothing more.

## 10. Independence vs. sum — the direct measurement
Take two individually-feasible b-tree bits (x_1530, x_1603) and their exact single-bit pin
solutions s1, s2 (each verified: applying it alone leaves ZERO atoms outside the selector core).

* **Shared variables:** s1 and s2 move 3 variables in common — `x_14853, x_31339, x_6083` —
  and they require **conflicting values on all three** (0 of 3 agree).
  Across all 50 feasible b-bits the 50 required values of `x_14853` are pairwise distinct,
  and distinct even mod p.
* **Not additive:** with both bits on and both solutions applied, the product flag
  `x_24195 = x_33953 * x_1250` flips from 0 to 1, and the two quantities
  `U = x_29210`, `V = x_8736` that the new rows constrain are **not** the sum of their
  singleton values, neither over Z nor mod p.

So the pin conditions are (a) NOT independent across bits — they collide on three shared
accumulators — and (b) NOT a function of the sum either; the coupling enters through the
*product* gates of the OR-tree.  Single-bit solutions therefore cannot be recombined
linearly, and each subset poses a fresh rank-2 congruence system in (U,V).

## 11. Single bits cannot close the selector either
For the branch that switches only the b-side on (`x_14853 = x_13682`), the closure was solved
at supports up to 3,811 variables / 2,727 rows (`runs/big01.log`) — **core infeasible at every
support**.  The blockers are, in raw form:
  a20215: `x_24530 = x_5647 * x_24908` with `x_5647 = 1`, `x_24530 = C1` -> demands `x_24908 = C1`;
  a28647: `x_36433 = x_26386*x_6083 + x_27475*x_33708` with `x_36433 = C2`,
with C1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
and C2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002.
Both are "produce this specific 296-bit constant on an accumulator" conditions.

## 12. Deliverables from the scan
`bitsol_<bit>_<score>.json` — one file per bit whose pin system solves *and* verifies with no
residual atoms outside the selector core.  Spot-checked with `../checker.py`:
`bitsol_10428_39017.json` -> **39017/39033** exactly, as the filename says.
