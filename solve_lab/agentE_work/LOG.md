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
(p = secp256k1 prime).  Fixing one pin relocates the residual one link down a chain
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
