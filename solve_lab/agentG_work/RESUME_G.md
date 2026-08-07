# Agent G — RESUME

## Best verified score
**39,026 / 39,033** — the inherited `solve_lab/best/new_instance_partial_39026.json`,
re-verified by me:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
I did not beat it, so nothing was written to `best/`. All my files are in
`solve_lab/agentG_work/`; no shared file was modified.

## The machine: exact symbolic forward evaluation over F_p
`gsym.py` (dense monomials) / `gsym2.py` (sparse monomials, scales to thousands of
symbols). Every gate output coefficient is ±1, so forward evaluation divides by nothing
and free-inputs -> every atom is an honest polynomial over Z, hence over F_p. Pick a set
of free inputs as indeterminates, evaluate every gate in topological order symbolically,
then every check atom.
* From `s10/AG_39013.json` with its booleans fixed: closed non-boolean symbol set = **112
  symbols**; of the **10,792 check atoms only 57 are non-constant** — 50 linear, 2
  quadratic, 5 cubic, **196 monomials total**. Validated at random points of F_p^112:
  **0 mismatches on all 10,792 checks**.
* Maximal model: ALL **6,117** non-boolean free inputs symbolic — 0.5 s pass, 0 gates
  skipped; 2,029 non-constant checks; sparse F_p elimination gives rank **1470**, 4,647
  free parameters, consistent; substitution makes every nonlinear check a CONSTANT, five
  nonzero — the same five values. The residual is exactly, algebraically pinned over the
  whole continuous freedom of the instance.
* EQUATION level (strictly more permissive than requiring each atom to vanish):
  **6,774 non-trivial equations = 6,613 linear (rank 1470, consistent) + 161 nonlinear**;
  forcing every linear equation leaves exactly **20** nonzero equations = AG_39013's 20
  failing lines, exactly.
* Validated at the deliverable: detaching the five gate outputs 7068/28730/29854/31864/642
  the model reproduces **0 atom mismatches over all 42,267 atoms and exactly 7 nonzero
  equations** at the deliverable's own point.

## The residual, in the instance's own variables
With K = 97553848499418123410591666447050222001188385549510401465815187079080512838891,

    A = (x22649 - x14853)^2 * (x22162 + x22649 + x14853 + K) - (x31339 - x16742)^2
    B = x14853*x30213 - x22649*x30213 + x14853*x16742 - x22649*x31339
        + x22162*x31339 - x16742*x22162
    a19297 = 8646263*A + 1073965*B ,  a19299 = 10159099*A + 6926539*B ,
    a30984, a36185, a40812 = three further members of the same rank-2 pencil.

The 50 linear checks pin all six of those variables to constants built from four literal
constants of the file, and there A != 0 and B != 0.

## Priority-1 answer: the boolean map is NOT affine on (A,B)
The 256 boolean free inputs that carry large load pins do not act affinely on (A,B) over
F_p. Setting one pins a wire to a specific ~256-bit literal AND re-routes which wire feeds
one of the six residual variables; two booleans acting on the same residual variable give
a value that is neither one's, nor the base, nor their sum. Since A is cubic and B
quadratic in those variables, bits -> (A,B) is a degree-3 polynomial map, not an affine
form. **There are no deltas (dA_i, dB_i) to sum, so the two-dimensional modular
subset-sum / LLL route does not apply.** Measured, not assumed.

## Minimum-weight coset decoding at equation level (the current task)
Model: 6,614 linear + 161 nonlinear equations in 6,122 unknowns.
* At the deliverable **all 7 violated equations are linear, zero nonlinear are violated**.
* Only **1,475 of the 6,122 unknowns occur in any linear equation and the linear rank is
  exactly 1,475** — full column rank, so the linear system pins every occurring unknown
  uniquely (point x0) and any departure costs equations.
* Cheapest unknowns: x22162 occurs in 2 linear equations {133,8073}; x30213 in 3; x9118
  and x29854 in 7; x8731 and x31864 in 9; x642 in 10. The deliverable's departure moves 15
  unknowns whose footprints union to 65 equations, 58 of which cancel.
* Departure on {x22162,x30213} alone: best point leaves **16 failing** (39,017).
* On the deliverable's 15-unknown support: 65 affine + 21 higher-degree equations; the
  affine rows have rank 11, giving a **4-dimensional cost-free departure space**
  (x1329, x9413, x10903, x17325) — on which all 20 cubics remain nonzero constants.
* Region closure: exactly **13** unknowns have their whole linear footprint inside the
  region; with x22162 and x30213 that is a closed **17-unknown** support, on which the 68
  affine rows collapse to **19 distinct directions** (multiplicities 1x15, 11, 13, 13, 16)
  of rank 13. Every violated-set of size <= 6 is therefore a subset of the 15
  multiplicity-1 directions.
* **EXHAUSTIVE result (`g66_exhaust17.py`, budgets 1..6, 4,880 admissible relaxations
  tested): NO relaxation of <= 6 equations lets the cubics be zeroed — every one leaves at
  least one cubic pinned to a nonzero constant.** So inside the closed region 7 is exactly
  optimal, established by exhaustion over the exact polynomial system rather than by
  search or by a tangent-space rank count.

## Re-enter
```
cd /home/user/integer_solver/solve_lab/agentG_work
python3 g11_bigsys.py                    # 112-symbol system + random-point validation
python3 g14_print.py                     # print the 57 polynomials explicitly
python3 g23_allsym.py ; python3 g24_bigsolve.py base
python3 g35_eqsolve.py -                 # equation-level exact solve
python3 g54_cosetsetup.py                # -> coset_model.pkl (the decoding instance)
python3 g56_colweight.py                 # per-unknown equation costs, x0, deliverable departure
python3 g66_exhaust17.py $(cat extsup.txt) 6    # exhaustive coset decoding, budget 6
```
NOTE `s9/eff/lib.py` does `os.chdir(solve_lab/s9)`; write outputs to absolute paths.

## Highest-value next experiment
The exhaustive optimality above is for the closed 17-unknown region. The one gap left is
whether a departure that also moves unknowns OUTSIDE that region (paying part of their
larger footprints but buying cancellation) can reach 6. Concretely: extend the support one
tier at a time — the 85 unknowns with 10 linear equations, then the 153 with 11 — and rerun
`g66_exhaust17.py` at budget 6. The direction-collapse (68 rows -> 19 directions) is what
makes the enumeration cheap, so the practical limit is how fast the number of
multiplicity-<=6 directions grows as the support widens.
