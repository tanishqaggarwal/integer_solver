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
