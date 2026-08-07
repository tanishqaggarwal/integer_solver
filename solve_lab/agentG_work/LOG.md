# Agent G LOG

## t0 — orientation
- Read PROMPT.txt, RESUME.md (1129 lines), STATE.json.
- Verified `best/new_instance_partial_39026.json` -> 39026/39033, 7 failing. CONFIRMED.
- Env: python-flint 0.9.0, sympy, z3, pysat, numpy. No scipy/sage/msolve/networkx/gmpy2.

## Exp G01 — coordinate/curve probe (`g01_probe.py`)
Prior claim (Part XXVI) says the core is secp256k1 point addition with
x1=x12186,y1=x16742,x2=x14853,y2=x24908,x3=x22162,y3=x30213.
Measured at 4 states:
- 39026 deliverable: x1==x2, y1==y2 exactly (mod p) -> the "addition" is degenerate there.
- AG_39013 / EC_39014 / PF_39015: P1 != P2.
- **None of P1,P2,P3 satisfies y^2 = x^3 + 7 (mod p) at any state.**
=> The EC reading is a formal algebraic identity, not literal secp256k1 points.
   Do NOT expect curve-theoretic structure (order n, discrete logs) to apply.

## Exp G02–G12 — EXACT SYMBOLIC REDUCTION OF THE WHOLE INSTANCE (the main result so far)
`gsym.py` performs **exact symbolic forward evaluation over F_p**: choose a set S of free
inputs as indeterminates, evaluate every gate in topological order as a genuine
multivariate polynomial (legal because every gate output coefficient is ±1), then evaluate
every CHECK atom symbolically.

Base state = `s10/AG_39013.json` (advice fixed point, score 39013), booleans fixed there.

* Support closure over non-boolean free inputs (`g10_closure.py`) -> **112 symbols**, closed.
* Symbolic pass (`g11_bigsys.py`): of the **10,792 check atoms, only 57 are non-constant**;
  the other 10,735 are the ZERO polynomial. Degrees: 50 linear, 2 quadratic, 5 cubic.
  **196 monomials in total.**
* VALIDATED EXACTLY: at 2 independent random points in F_p^112 the symbolic prediction
  matches the true forward-evaluated value of all 10,792 checks — 0 mismatches
  (`g07_validate.py`, `g11_bigsys.py`).
* Linear part: 50 equations, **rank 38, consistent, 74 free parameters**.
  Substituting the general linear solution into the 7 nonlinear checks
  (`g12_solve112.py`) makes every one of them a **CONSTANT**: two are 0, five are nonzero.
  => in this boolean frame the residual is **EXACTLY, ALGEBRAICALLY infeasible mod p** —
  not a tangent-space statement.

### What the 5 residual checks actually are (`g14_print.py`)
All five are the same rank-2 pencil in
  x1=x22649, y1=x16742, x2=x14853, y2=x31339, x3=x22162, y3=x30213:
  A' = (x2-x1)^2 * (x3 + x1 + x2 + K) - (y2-y1)^2
  B  = (y3+y1)*(x2-x1) - (x1-x3)*(y2-y1)
  a19297 = 8646263*A' + 1073965*B, a19299 = 10159099*A' + 6926539*B, etc.
So the entire instance reduces to **A' = 0 and B = 0 mod p**. (Confirms Part XXVI/XXVII's
identities, and my derivation is independent and exact.)
The other 50 linear checks pin all six coordinates to constants through 4 literal pins:
  a1618 -> x3 ; a688 -> y3 ; a31670(x24601) -> x22152 -> a2423 -> x1 ;
  a3576(x2081) -> x6418 -> a29539 -> x2 ; a31672(x24601) -> x33462 -> x8778 -> y1 ;
  a3578(x2081) -> x12553 -> x24548 -> x14623 -> y2.

### Branch table, computed EXACTLY (`g18_combo.py`), not by local repair
| flip | ninc(linear) | nonzero-const checks | nonlinear residual |
|---|---|---|---|
| (base) x2081=1,x24601=1 | 0 | 0 | **5 (A'!=0,B!=0)** |
| x2081=0 | 5 | 0 | 0 |
| x24601=0 | 3 | 0 | 0 |
| both 0 | 0 | **6** | 0 |
In the (0,0) branch the six unreachable checks are a688, a1618 (the x3/y3 pins, now
demanding 8863713*x18956 = C with x18956 stuck at 0 because the selector x15298 = 0),
a23000 (the forced OR gate x9274 = 1), a39067, a40608, a41211. This is the covering
design, now visible as exact algebra rather than as a search outcome.

## Exp G23–G38 — the maximal sound model, and what it settles
* `gsym2.py` = sparse-monomial version of the symbolic evaluator. With **ALL 6,117
  non-boolean free inputs** as symbols the symbolic pass takes **0.5 s** and produces
  **2,029 non-constant check polynomials** (1,883 linear, 141 quadratic, 5 cubic),
  **0 gates skipped**, 0 unreachable nonzero constants.
* Sparse F_p elimination (`gsolve.py`): rank **1470**, 4,647 free parameters, consistent.
  Substituting into the nonlinear part leaves the SAME five nonzero constants as the
  112-symbol closure. **So the mod-p infeasibility of the AG_39013 boolean frame holds
  over the entire continuous freedom of the instance — no support approximation, no
  dead-monomial blindness.**
* EQUATION-level version (`g35_eqsolve.py`, strictly weaker than atom-level): 6,774
  non-trivial equations, 6,613 linear (rank 1470, consistent) + 161 nonlinear;
  forcing every linear equation leaves **exactly 20 nonzero equations** = precisely
  AG_39013's 20 failing. Exact, not searched.
* MODEL VALIDATION AT THE DELIVERABLE (`g38_delivcheck.py`): detaching the five gate
  outputs 7068/28730/29854/31864/642, my symbolic model evaluated at the 39,026
  deliverable's own point reproduces **0 atom mismatches over all 42,267 atoms and
  exactly 7 nonzero equations mod p**. The deliverable lies inside my parameter space;
  it beats the "all linear equations forced" point (20) by deliberately violating
  linear equations, i.e. it is a minimum-weight coset leader, not a linear optimum.

### The MUX, mapped exactly (`g28_showcheck.py`, `g29_frame.py`, `g32_framesolve.py`)
`a1618` and `a688` are a 2-bit multiplexer on which coordinate pair they pin:
| (x2081,x24601) | a1618 pins | a688 pins | core a19297 |
|---|---|---|---|
| (1,1) base | x22162 = x3 | x30213 = y3 | LIVE: 5 residual constants (A!=0,B!=0) |
| (0,1) | x22649 = x1 | x16742 = y1 | dead; 5 inconsistent linear certs |
| (1,0) | x14853 = x2 | x31339 = y2 | dead; 4 inconsistent linear certs |
| (0,0) | CONSTANT (unsatisfiable) | CONSTANT | dead; 6 unreachable checks |
Every branch clashes with the coordinate's own pin chain. With x4287 or x13195 on, the
residual becomes POLYNOMIAL in (x14853,x31339) = (x2,y2) — but the branch obligation
checks (a19088/a22233/a22235 for x4287; a7932/a7934/a7936/a41512 for x13195) pin those
two variables to a unique point at which A,B != 0. `g31_solvres.py` reproduces the
session-9 "branch obligation" residue 33371159155735472537534252650716501592825364489306217536352743247010353604716
independently and exactly.
Frame 4287+13195 gives 4 unknowns (x8731,x9118,x14853,x31339) and 13 polynomials;
`g33_solve2.py` eliminates all four and leaves 4 nonzero residuals. Not a solution.
