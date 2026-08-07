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

## Exp G39–G40 — the message bits REWIRE the coordinates (partial; scans stopped)
`g39_pinscan.py`: flipping one message bit changes WHICH WIRE feeds a coordinate.
e.g. bit x47 moves y1's chain from x8778 to x8060 and x1 from x22649 to x28548, both of
which are pinned to 0 by a32226 / a20109. So the message selects the six coordinates from
a MUX tree of pinned wires — and selection and pinning are the SAME bit, which is the
trapdoor: a wire cannot be selected and left free at the same time.
`g40_affine.py` (written, not completed): tests whether the bits act affinely on (A,B).
That measurement is the input to the LLL / low-density-subset-sum attack in RESUME_G.

## STOP checkpoint (coordinator)
`g27_bigscan.py` (full exact per-boolean reduce) reached ~250/1156; `g39_pinscan.py`
~300/1156; both terminated. Partial pickles: `bigscan_all.pkl`, `pinscan.pkl`.
Best verified score remains the inherited 39,026 — I did not beat it.
All my artifacts are under solve_lab/agentG_work/; no shared file was modified.

## Exp G41 — MY REFUTATION (i) IS WRONG. RETRACTED.
I tested the curve in the wrong (short) form. Re-run in general Weierstrass form,
deriving every constant from my own A,B identities and nothing else (`g41_curve.py`):
* On `y^2 = x^3 + a2 x^2 + a4 x + a6` the chord-and-tangent addition is exactly
  `A = 0` with **K = a2**, and `B = 0`. So the constant K in my A identity IS a2.
  a2 = K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
* Solving for a4, a6 from P1, P2:
  a4 = 114170008767671698752186727197936107864370654164657728518655355473804451402762
  a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490
* **a4 == a2^2/3 exactly, so j = 0.** After x -> x - a2/3: A_short = 0 (exactly),
  B_short = 64019533680030876408443198762210829058751700634554282185987325820393598524794.
  B_short/7 is a **6th power** mod p -> the curve is ISOMORPHIC to secp256k1 over F_p.
* P1, P2 **and** P3 all lie on it; `[n_secp]P1 = [n_secp]P2 = O`.
* Isomorphism: `x_sec = (x + a2/3)/u^2`, `y_sec = y/u^3`,
  u = 4210889811980686189396764679825672592540066047176031544704936155054310740018.
Agent C's finding is confirmed independently. My earlier negative came from assuming a2=0.

## Exp G42–G51 — WHAT THE INSTANCE ACTUALLY IS
* `g42_secp.py`: the three coordinates map to genuine secp256k1 points of order n.
  P1+P2 != P3, so the residual `A=B=0` is exactly the assertion **P3 = P1 + P2**.
* `g46_table.py` / `g45_pointscan.py`: flipping a message bit FREES a coordinate pair and
  lights load-pin obligations that re-pin it — and the re-pinned value is again a genuine
  curve point. Each bit therefore selects a table point.
* `g47_pairs.py`: two bits acting on the SAME coordinate give a value that is neither
  bit's point, nor the base, nor their group sum, and is OFF the curve — same-coordinate
  bits conflict. Cross-coordinate pairs are group-affine in D = P3-(P1+P2), trivially.
* `g49_loads.py`: **exactly 256** boolean free inputs carry huge load pins (886 pins);
  602 of the 886 constants are valid x-coordinates vs ~443 at random.
* `g50_points.py`: pairing each bit's load constants into curve points and testing the
  doubling map gives **255 doubling hits over 256 points**.
* `g51_chain.py`: the 256 points form ONE chain with a unique root; ordering it and
  recomputing gives **P(bit_i) = [2^i]P0 for all 256 positions, 0 exceptions**.
  P0 is a secp256k1 point of order n (not G, not a small multiple of G).
* `g53_export.py` (re-verified, `secp_identification.json`):
  **P1 = [2^72]P0 (bit x24601), P2 = [2^235]P0 (bit x2081)** — the two bits set in the
  base state — and P3 is the pinned target, on the curve, of order n.

### THEREFORE
The message is a 256-bit scalar `k = sum_i b_i 2^i`; the circuit is a double-and-add of
`[k]P0` and the instance asserts `[k]P0 = P3`. With only two bits on, the double-and-add
degenerates to a single addition — which is exactly why my exact reduction found ONE
point addition and 57 non-constant checks. **A full solve is precisely
`k = log_{P0}(P3)` on secp256k1.**
This settles the coordinator's dichotomy on the NEGATIVE side: the bits do NOT act
affinely on (A,B) — they act by GROUP DOUBLING. `sum b_i [2^i]P0 = P3` is the discrete
logarithm itself, not a low-density knapsack (density 1, and the map is the group
exponential, not an F_p-linear form), so no LLL / subset-sum attack applies.

## Exp G52 — the one cheap shot at the trapdoor
`g52_lowweight.py`: meet-in-the-middle over subsets of the 256 chain points.
Hamming weight of k **<= 4 is ruled out** (weight 1,2 exhaustive; 3,4 by MITM);
the weight<=6 pass was still running at cutoff. `k = 2^72 + 2^235` (the current state)
is confirmed NOT a solution. Nothing below weight 5 works, so k is a generic 256-bit
scalar and Pollard rho (~2^128) is the only generic route. I stop here honestly:
the instance is satisfiable by construction (the setter knows k) but finding k is ECDLP.
