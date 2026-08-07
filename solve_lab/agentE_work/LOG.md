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

## 13. THE TRIPLE IS SOLVABLE — the recombination barrier falls (my §10 conclusion is REFUTED)

Target: the size>=2 obstruction from §9.4, at the repaired pair state (x_1530, x_1603),
`triple_state_seed.json`, whose only non-selector bad atoms are {722, 724, 726}.

### 13.1 Reduction of the three conditions to two
`p | U` and `p | V` and `5002401*U + 15322661*V = 0` collapse: if `p | U` then from the exact
row `5002401*U = -15322661*V` and `gcd(15322661, p) = 1` give `p | V` for free.  So the triple
is exactly **{ p | U ,  5002401*U + 15322661*V = 0 over Z }**, then two quotient handles.

### 13.2 The handles the linear model never proposed
Of the 23 free variables that move (U,V), **14 are boolean-constrained** (an atom
`b*(b-1)`, `b*b-b` or `2*b*(1-b)` pins them to {0,1}) and **9 are unconstrained integers**
(x_3401, x_4012, x_25710, x_28954, x_30468, x_33169, x_33177, x_34801, x_37856).
Of those 9, `(U,V)` is **exactly affine** in **x_30468 and x_33169** (checked at +1, +2, +7);
`dU/dx_33169 = 0`, so the pair is triangular.  Earlier searches never proposed them because
they sat in rows the closure had dropped as nonlinear.

### 13.3 Exact two-unknown congruence solve (`triple4.py`)
    d1 = 107833348867425503451660448077443469436092985532849038009363238620095543518536 + k*p
    d2 = -(N0 + M*d1)/(15322661*D)     with M = 5002401*A + 15322661*B
and the second congruence has modulus 1, so **every** k works — a one-parameter family.
Result, exactly re-evaluated: `p|U = True`, `p|V = True`, `5002401*U + 15322661*V = 0 = True`;
atom 726 vanishes, and setting the two free handles `x_34496 = x_6635/p`,
`x_3193 = 7952523*x_36280/p` kills 722 and 724 as well.  **The triple closes.**

### 13.4 Simultaneous affine solve (`triple7.py`, `triple8.py`)
Rebuilt as one system: 56 knobs on which the *whole residual vector* is affine (verified at
+1,+2,+7) x 46 atoms.  Full system: infeasible.  Greedy maximal solvable subset:
**44 of 46 rows solve simultaneously**, and the only two blocking rows are

    a20215 : x_24530 - x_5647*x_24908
    a28647 : x_36433 - (x_36990 + x_19239)

Applying the 44-row solution and re-evaluating exactly gives **39,005/39,033 with just TWO
nonzero atoms in the entire instance** (`triple8_39005.json`, verified below).

**So the triple is NOT the barrier.**  What binds is a20215 and a28647 — the same two rows
that blocked the single-b-bit branch in §11.  Both are "an accumulator must equal a specific
296-bit constant" conditions:
  a20215 -> `x_24908 = C1`,  a28647 -> `x_26386*x_6083 + x_27475*x_33708 = C2`.

### 13.5 Verification note
The solution's values reach 4,430 decimal digits, above Python's default 4,300-digit
string-conversion cap, so `checker.py` cannot *parse* the file with default settings.
`verifyE.py` raises only that cap and then calls `checker.load_equations`,
`checker.load_assignment` and `checker.evaluate_all` unmodified:
    python3 verifyE.py triple8_39005.json  ->  satisfied 39005/39033 (28 failing)

## 14. a20215 / a28647 — the affine solve, and the exact reason it stops
Knob set built as instructed: every free variable in the cones of x_24908 and of
(x_6083, x_33708), plus the cones of the four cluster atoms, **with the non-boolean integer
handles explicitly included** (the inclusion that unlocked the triple).

### 14.1 The cluster is a closed 5-row system
From `triple8_seed.json` (bad = {20215, 28647}) six integer knobs have their entire
disturbance set inside {7389, 10187, 20212, 20215, 28647}; adding the two non-boolean movers
of a20212 gives eight.  Exact rows (`close5.py`, `runs/close7.log`):

    a7389 :   d_6083  - p*d_26489                       = 0
    a10187:   d_31339 - c1*d_37012                      = 0
    a20212:   p*d_11436 + c2*d_14393 - d_14853          = 0
    a20215: - p*d_22820 - d_31339                       = R1
    a28647: - d_6083 + d_14853                          = R2

R1 = resid(a20215), R2 = resid(a28647).

### 14.2 Why it is infeasible — measured, not assumed
    R1 mod p = 22981624690591324143788809642515852940280603493270692712106986169263210356252
    R2 mod p = 44159679639019146557987083382852396884224992023970032213706899677695745279353
    c1 = d(a10187)/dx_37012 : 279 bits, **c1 = 0 (mod p)**
    c2 = d(a20212)/dx_14393 : 280 bits, **c2 = 0 (mod p)**
So a10187 forces `d_31339 = 0 (mod p)` while a20215 wants `d_31339 = -R1 (mod p)`;
a20212 forces `d_14853 = 0 (mod p)` while a28647+a7389 want `d_14853 = R2 (mod p)`.
**Every knob that reaches a10187 or a20212 enters with a coefficient divisible by p**, so both
chains are mod-p trivial and the two non-zero residues R1, R2 cannot be absorbed.
This is the same shape as "a row pins a knob to zero" — an independent second instance of it.

### 14.3 What would break it
A knob reaching a10187 or a20212 with coefficient **coprime to p**.  Searched: the full cones
of x_24908, x_6083, x_33708, x_24530, x_36433, x_36990, x_19239, x_26386, x_27475, x_5647 and
of the four cluster atoms — 326 candidates, 248 cluster-affine, only 8 with a contained
disturbance set, and among the a20212 movers only 3 are non-boolean (x_11436 with coefficient
exactly p, x_14393 with c2 = 0 mod p, x_14853 with coefficient -1 but it is the unknown itself).
The 181 remaining movers of a20212 are all boolean selector bits.

## 15. The 181 boolean movers, as 0/1 DECISIONS — §14's claim is FALSE as stated
Measured by actual re-propagation at 0 and 1 (not derivatives), from `triple8_seed.json`
(`boolknob.py`, `runs/boolknob.log`).  Census of the exact residual delta
`Delta_b = resid(b=1) - resid(b=0)` over the 256 boolean movers in the two cones:

| row | delta = 0 | delta = 0 (mod p) | **delta NONZERO mod p** |
|---|---|---|---|
| a10187 | 233 | 0 | **23** |
| a20212 | 78  | 0 | **178** |

**So a knob reaching a10187 and a20212 with content coprime to p DOES exist — 201 of them.**
LOG §14's "every knob reaching a10187 or a20212 enters with a coefficient divisible by p" was
measured on the AFFINE class only and generalised past the evidence.  Corrected here.

### 15.1 But counting cannot close it
The nonzero deltas are nearly all the *same* residue:
  a20212: **one** distinct residue, 36200939269128454586076546451607958467047992891178506183612554289882454126226, x178
  a10187: **three** residues, x21 / x1 / x1
Required multiplicities to absorb R1, R2 (from the §14 rows) are full 256-bit numbers:
  k = 101108319720394122322776115727804160327749809720865796405265564477795935003136 (a10187)
  m = 106495465405897155704077983448098605801616252488421425517764957547039336632879 (a20212)
so no subset of <=23 (resp. <=178) equal-residue flips supplies the needed shift by counting.
Composition experiments (`boolsolve.py`, `boolfix.py`): flipping a bit brings in that bit's own
pin atoms, and the 8-knob cluster solve then fails on the *pin* rows, not the cluster rows;
flip + full iterated pin repair reaches at best the same 39,005.

### 15.2 Context dependence — the same signature a third time
In the affine picture at this state every knob on a10187/a20212 was mod-p trivial.  With two
selectors already on, **178 booleans become mod-p nontrivial**.  The coefficient structure
therefore genuinely depends on which selectors are on; the mod-p content of these rows is not
a fixed property of the row.
