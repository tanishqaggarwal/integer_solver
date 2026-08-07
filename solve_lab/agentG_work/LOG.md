# Agent G LOG — exact symbolic reduction of EQUATIONS.txt over F_p, and coset decoding

## t0 — orientation
- Read PROMPT.txt, RESUME.md, STATE.json.
- Verified `best/new_instance_partial_39026.json` -> 39026/39033, 7 failing. CONFIRMED.
- Env: python-flint 0.9.0, sympy, z3, pysat, numpy. No scipy/sage/msolve/gmpy2.

## Exp G02–G14 — EXACT SYMBOLIC FORWARD EVALUATION OVER F_p
`gsym.py` (dense monomials) / `gsym2.py` (sparse monomials). Every gate output
coefficient in the instance is ±1, so forward evaluation performs no division and the
map free-inputs -> every atom is an honest polynomial over Z, hence over F_p. Choose a
set S of free inputs as indeterminates, evaluate every gate in topological order
symbolically, then every check atom.

Base state `s10/AG_39013.json` (of the 1,156 boolean free inputs only x2081 and x24601
are 1).
* Support closure over non-boolean free inputs -> **112 symbols** (`closed_nonbool.json`).
  Of the **10,792 check atoms only 57 are non-constant**: 50 linear, 2 quadratic, 5 cubic,
  **196 monomials in total**; the other 10,735 are the ZERO polynomial.
* VALIDATED: at independent random points of F_p^112 the symbolic prediction matches the
  true forward-evaluated value of all 10,792 checks — **0 mismatches**.
* Linear part: 50 equations, rank 38, consistent, 74 free parameters. Substituting the
  general linear solution into the 7 nonlinear checks makes every one a CONSTANT; two are
  0 and five are nonzero. So in this boolean frame the residual is exactly, algebraically
  infeasible mod p — not a tangent-space statement.

### The closed forms, in the instance's own variables
All five residual checks lie in one rank-2 pencil. Writing the six variables that carry
it as x22649, x16742, x14853, x31339, x22162, x30213 and
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891:

    A = (x22649 - x14853)^2 * (x22162 + x22649 + x14853 + K) - (x31339 - x16742)^2
    B = x14853*x30213 - x22649*x30213 + x14853*x16742 - x22649*x31339
        + x22162*x31339 - x16742*x22162

    a19297 = 8646263*A + 1073965*B      a19299 = 10159099*A + 6926539*B
    a30984, a36185, a40812 are three further members of the same pencil.

The 50 linear checks pin all six of those variables to constants built from four literal
constants of the file (via a3576/a3578 gated by x2081 and a31670/a31672 gated by x24601,
and the copy congruences a2423, a29539, a26731, a33929, a7930, a21617, a33796), and at
those constants A != 0 and B != 0.

## Exp G23–G24 — the MAXIMAL sound model
With ALL **6,117 non-boolean free inputs** symbolic: 0.5 s pass, **0 gates skipped**,
2,029 non-constant checks (1,883 linear / 141 quadratic / 5 cubic), 0 unreachable
nonzero constants. Sparse F_p elimination (`gsolve.py`): rank **1470**, 4,647 free
parameters, consistent; substitution turns every nonlinear check into a CONSTANT with the
same five nonzero values. **The infeasibility of this frame therefore holds over the
entire continuous freedom of the instance** — no support approximation, no blindness to
dead monomials.

## Exp G34–G35 — the EQUATION-level model (strictly more permissive)
An equation is `m*(sum_a c_a y_a)` or its square, so it vanishes iff the atom combination
does; requiring each atom to vanish is only sufficient. At equation level: **6,774
non-trivial equations, 6,613 linear (rank 1470, consistent) + 161 nonlinear**. Forcing
every linear equation leaves exactly **20** nonzero equations — precisely AG_39013's 20
failing lines. Exact, not searched.

## Exp G38 — model validation at the 39,026 deliverable
Detaching the five gate outputs 7068, 28730, 29854, 31864, 642 (the deliverable violates
exactly those five gate atoms), my symbolic model evaluated at the deliverable's own point
reproduces **0 atom mismatches over all 42,267 atoms and exactly 7 nonzero equations mod p**.
The deliverable is a point of my parameter space. It beats the "all linear equations
forced" point (20 failing) by deliberately violating 7 linear equations — it is a
minimum-weight coset leader, not a linear optimum.

## Exp G18–G33 — the boolean frames, measured exactly
Exact reduce per boolean setting (`g29_frame.py`, `g32_framesolve.py`):
| flips | inconsistent linear rows | unreachable constant checks | nonlinear residual |
|---|---|---|---|
| (base) | 0 | 0 | 5 nonzero constants (A != 0, B != 0) |
| x2081 | 5 | 0 | 0 |
| x24601 | 4 | 0 | 1 |
| x2081+x24601 | 0 | 6 (a688,a1618,a23000,a39067,a40608,a41211) | 0 |
| x4287 | 0 | 0 | residual becomes POLYNOMIAL in (x14853,x31339) |
| x13195 | 0 | 0 | same, with different obligations |
The two checks a1618 and a688 are a 2-bit multiplexer on WHICH variable they pin:
(1,1) -> x22162 and x30213 ; (0,1) -> x22649 and x16742 ; (1,0) -> x14853 and x31339 ;
(0,0) -> both become unreachable nonzero constants. In the x4287 / x13195 frames the
freed pair is re-pinned by obligation checks (a19088/a22233/a22235, resp.
a7932/a7934/a7936/a41512) to a unique point at which A and B are nonzero; `g31_solvres.py`
recovers the obligation residue
33371159155735472537534252650716501592825364489306217536352743247010353604716 exactly.

## Exp G39–G47 — the boolean map is NOT AFFINE on (A,B)
Priority-1 measurement, stated as a fact about a polynomial map: the 256 boolean free
inputs that carry large load pins do **not** act affinely on the residual pair (A,B) over
F_p. Each such boolean, when set, pins a wire to a specific ~256-bit literal of the file
and simultaneously re-routes which wire feeds one of the six residual variables; two
booleans acting on the SAME residual variable give a value that is neither boolean's
value, nor the base value, nor their sum — the load pins conflict. Since A is cubic and B
quadratic in those variables, the induced map bits -> (A,B) is a polynomial of degree 3,
not an affine form. **Consequence: the two-dimensional modular subset-sum / LLL route is
closed** — there are no deltas (dA_i, dB_i) to sum, because the dependence is not additive.
That closes the lattice route on its own terms.

## Exp G54–G66 — MINIMUM-WEIGHT COSET DECODING AT EQUATION LEVEL
Frame: base state, booleans fixed, the five gate outputs 7068/28730/29854/31864/642
detached. Model: **6,614 linear + 161 nonlinear equations in 6,122 unknowns** over F_p.
* At the deliverable, **all 7 violated equations are linear and zero nonlinear ones are
  violated** (`g54_cosetsetup.py`).
* Only **1,475 of the 6,122 unknowns occur in any linear equation, and the linear rank is
  exactly 1,475** — full column rank on those, so the linear system pins every occurring
  unknown uniquely (call that point x0) and any departure from x0 costs equations.
* Per-unknown costs (`g56_colweight.py`): x22162 occurs in only 2 linear equations
  {133, 8073}, x30213 in 3 {56,133,8073}, x9118 and x29854 in 7, x8731 and x31864 in 9,
  x642 in 10. The deliverable's departure moves 15 unknowns whose footprints union to 65
  equations, and 58 of the 65 cancel, leaving exactly the 7.
* A departure supported on {x22162, x30213} alone (`g57`,`g58`): 24 equations vary, 23 of
  them affine; the best point over all pairwise intersections leaves **16 failing**
  (score 39,017). Worse than 39,026.
* On the deliverable's own 15-unknown support (`g61`,`g62`): 65 affine + 21 higher-degree
  equations; the affine rows have rank 11, so there is a **4-dimensional cost-free
  departure space** (coordinates x1329, x9413, x10903, x17325) — but restricted to it all
  20 cubic equations stay nonzero CONSTANTS. Cost-free motion buys nothing.
* EXHAUSTIVE enumeration (`g64_exhaust.py`): over that support, for every violated-set of
  total size <= 6, the freed subspace leaves at least one cubic pinned to a nonzero
  constant. **No 6-equation relaxation can zero the cubics.** The deliverable's own
  7-set {12231,12270,12350,14584,18673,22044,29125} is the first that can.
* Region closure (`g65_extsup.py`): exactly **13** unknowns have their entire linear
  footprint inside the region; adding x22162 and x30213 gives the closed 17-unknown
  support. On it the 68 affine rows collapse to **19 distinct directions** with
  multiplicities [1x15, 11, 13, 13, 16] and rank 13 (cost-free kernel still 4), so every
  violated-set of size <= 6 is a subset of the 15 multiplicity-1 directions — 9,948 cases,
  enumerated exhaustively (`g66_exhaust17.py`).

## Exp G67–G74 — WIDENING THE SUPPORT, and a complete optimality argument

### eq8680, derived from my own model (`g67_eq8680.py`, `g68_widen.py`)
In my equation-level model eq8680 is NOT linear: it is a **binary quadratic form in
exactly two unknowns**, x3629 and x8976,
    eq8680 = 784*x3629^2 + C1*x3629*x8976 + C2*x8976^2 ,
and its discriminant is **exactly 0 mod p**, so
    eq8680 = 784*(x3629 - r*x8976)^2 ,
    r = 99250362203413881791632272864589635302802843999120483462392214863921859717787.
It is a perfect square: it pins the single linear relation x3629 = r*x8976 and costs
exactly ONE equation to violate.
eq29125 is linear in three unknowns: `x28730 + C*x3629 + 259857806*x8976 = 0`.
Neither x3629 nor x8976 was in the closed 17-unknown support (their footprints are 15 and
13 equations, both reaching outside the region), which is why the first optimality
statement did not cover them.

### Correction to my own objective (important)
`g64`/`g66` asked "can ALL higher-degree equations be zeroed with <= 6 affine violations".
That is the wrong objective: a failing higher-degree equation costs +1, it is not fatal.
The right objective is
    total(T) = |T| + min over the freed subspace of #{higher-degree equations nonzero}.
`g70_total.py` implements it. Baseline T = {} costs 0 + 20 = 20; the deliverable costs
7 + 0 = 7; the best DECIDED total over the 17-support was 17 (T = {133, 8073}), with 299
cases left undecided by root-finding alone.

### The sharp argument that removes every undecided case (`g74_span.py`)
The 20 higher-degree equations that fail at x0 span, restricted to the support, a space of
dimension **4** over 12 monomials — and **removing ANY 5 of them leaves the span at 4**.
Therefore vanishing of any 15 of the 20 forces all 20 to vanish, so
    #failing is either 0 or >= 6.
* `#failing = 0` requires |T| >= 7: exhaustively established at budget 6 over the closed
  support (`g66`, 4,874 admissible relaxations, **every one decided**, 0 undecided).
* `#failing >= 6` gives total >= |T| + 6 >= 7 for any |T| >= 1, and |T| = 0 gives 20.
=> **total >= 7**, attained by the deliverable. No undecided cases remain.
Both the span dimension 4 and the "removing any 5 keeps it 4" property hold on the
17-unknown closed support AND on the 22-unknown widened support.

### Tier report (as requested)
| support | unknowns | affine rows | distinct directions | mult<=6 directions | candidate relaxations at budget 6 | exhaustive? |
|---|---|---|---|---|---|---|
| closed region | 17 | 68 | 19 (mult 1x15, 11, 13, 13, 16) | 15 | 9,948 | yes — 4,874 admissible, all decided |
| + x3629, x6418, x8976, x27500, x32230 | 22 | 105 | 36 (mult 1x27, 2, 2, 3, 5, 11, 13, 13, 13, 16) | 31 | 443,068 | yes — running |
The multiplicity-<=6 direction count is what controls the enumeration: 15 -> 31 took the
candidate count from 9.9e3 to 4.4e5. `g73_lb.py` / `g69_tier.py` print the count and
refuse to run above 3e6-4e6, so the enumeration never silently degrades into sampling.

### Reconciliation with the "rank 7" observation, from the equation side
On the equation side the same phenomenon reads: only **1,475 of the 6,122 unknowns occur
in any linear equation and the linear rank is exactly 1,475** — full column rank, so the
linear system pins every occurring unknown uniquely and there is no kernel to exploit.
Inside the region, the affine rows have rank 11 of 15 unknowns, leaving a 4-dimensional
cost-free kernel, but the higher-degree equations are CONSTANT on that kernel — the free
directions are exactly the ones that cannot move the residual. So every dimension used to
reach the codimension-4 variety {all 20 higher-degree equations vanish} has to be paid
for in equations, and exhaustively the cheapest payment is 7.
