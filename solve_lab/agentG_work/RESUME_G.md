# Agent G — RESUME (relaxation & rounding -> EXACT symbolic reduction over F_p)

## Best verified score
**39,026 / 39,033** — the inherited `solve_lab/best/new_instance_partial_39026.json`,
re-verified by me: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
I produced NO improvement; nothing of mine reaches 39,026, so nothing was added to best/.

## THE REUSABLE RESULT: exact symbolic forward evaluation over F_p
`agentG_work/gsym.py` (dense monomials) and `agentG_work/gsym2.py` (sparse monomials,
scales to thousands of symbols). Every gate output coefficient is +-1, so forward
evaluation divides by nothing and free-inputs -> every atom is an honest polynomial over
Z, hence over F_p. Pick a symbol set S of free inputs, evaluate every gate in topological
order symbolically, then every check atom.
Base state `s10/AG_39013.json` (booleans: only x2081=x24601=1, other 1154 are 0).
* Closed non-boolean symbol set = **112 symbols** (`closed_nonbool.json`).
  Of the **10,792 check atoms only 57 are non-constant**: 50 linear, 2 quadratic,
  5 cubic, **196 monomials total**; the other 10,735 are the ZERO polynomial.
  Validated at random points of F_p^112: **0 mismatches on all 10,792 checks**.
* Maximal model: ALL **6,117 non-boolean free inputs** as symbols -> 0.5 s pass,
  2,029 non-constant checks (1,883 lin / 141 quad / 5 cubic), **0 gates skipped**.
  Sparse F_p elimination: rank **1470**, 4,647 free params, consistent; substitution
  makes every nonlinear check a CONSTANT and five are nonzero — the same five values.
  => the mod-p infeasibility of this boolean frame holds over ALL continuous freedom
  (no support approximation, no dead-monomial blindness).
* Equation-level version (weaker, more permissive): 6,774 non-trivial equations,
  6,613 linear (rank 1470, consistent) + 161 nonlinear; forcing every linear equation
  leaves exactly **20** nonzero equations = AG_39013's 20 failing, exactly.
* The whole instance mod p is `A = 0 and B = 0`:
  `A = (x1-x2)^2 (x3+x1+x2+K) - (y2-y1)^2`
  `B = y3(x2-x1) + y1 x2 - x1 y2 + x3 y2 - x3 y1`
  x1=x22649, y1=x16742, x2=x14853, y2=x31339, x3=x22162, y3=x30213,
  K = 97553848499418123410591666447050222001188385549510401465815187079080512838891.
  a19297 = 8646263A+1073965B, a19299 = 10159099A+6926539B; a30984, a36185, a40812 are
  three more members of the same rank-2 pencil. All six coordinates are pinned to
  constants by the 50 linear checks, and at those constants A != 0, B != 0.
* MODEL VALIDATED AT THE DELIVERABLE (`g38_delivcheck.py`): detaching the five gate
  outputs 7068/28730/29854/31864/642, the model evaluated at the 39,026 point gives
  **0 atom mismatches over all 42,267 atoms and exactly 7 nonzero equations mod p**.

## Commands to re-enter
```
cd /home/user/integer_solver/solve_lab/agentG_work
python3 g11_bigsys.py                 # 112-symbol system + random-point validation
python3 g14_print.py                  # print all 57 polynomials explicitly
python3 g23_allsym.py [state] [flips] # maximal 6117-symbol model -> allsym_*.pkl
python3 g24_bigsolve.py base          # sparse linear solve + substitution
python3 g35_eqsolve.py -              # EQUATION-level exact solve
python3 g29_frame.py - 2081 24601 4287 13195 2081,24601
python3 g32_framesolve.py 2081,24601,4287,13195 4   # + solvability verdict per frame
python3 g38_delivcheck.py             # model validation at the 39,026 deliverable
```
NOTE: `s9/eff/lib.py` does `os.chdir(solve_lab/s9)`, so scripts writing relative paths
land in `s9/`. Use absolute paths or move the outputs back (I did).

## Single highest-value next experiment
Finish `g40_affine.py`: test whether the 1,156 boolean free inputs act AFFINELY on (A,B)
over F_p. If affine, collect all deltas (dA_i,dB_i) and solve the 2-dimensional modular
subset-sum `sum_i b_i d_i = -(A0,B0) mod p, b in {0,1}` by LLL — density 256/512 = 0.5
sits inside the low-density regime where lattice attacks succeed. Partial per-bit data:
`bigscan_all.pkl` (~250/1156 exact reduces), `pinscan.pkl` (~300/1156 rewiring signatures).
