# Agent G — RESUME (relaxation & rounding -> exact algebraic reduction)

## Best verified
**39,026 / 39,033** — the inherited `solve_lab/best/new_instance_partial_39026.json`
(re-verified by me: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`). No improvement of my own yet.

## Main result (new, exact, validated)
`agentG_work/gsym.py` does EXACT symbolic forward evaluation over F_p. From the state
`s10/AG_39013.json` with its booleans fixed, closing over non-boolean free inputs gives
**112 symbols**, and of the 10,792 check atoms only **57 are non-constant polynomials**
(50 linear, 2 quadratic, 5 cubic; 196 monomials). Validated at random points: 0 mismatches
on all 10,792 checks. Solving the linear part (rank 38, consistent, 74 free params) and
substituting turns every nonlinear check into a CONSTANT — 5 nonzero.
=> The whole instance mod p is exactly `A' = 0 and B = 0` with
   `A' = (x2-x1)^2 (x3+x1+x2+K) - (y2-y1)^2`, `B = (y3+y1)(x2-x1) - (x1-x3)(y2-y1)`,
   x1=x22649, y1=x16742, x2=x14853, y2=x31339, x3=x22162, y3=x30213,
   and all six coordinates are pinned to constants by the 50 linear checks.
Branch table (exact): (1,1) -> 5 nonlinear residuals; (0,1) -> 3 inconsistent linear rows;
(1,0) -> 5; (0,0) -> 0 inconsistent but 6 unreachable constant checks (a688,a1618,a23000,
a39067,a40608,a41211).

## Re-enter
```
cd /home/user/integer_solver/solve_lab/agentG_work
python3 g18_combo.py 2081 24601 2081,24601      # exact reduce for boolean flips
python3 g15_build.py <state.json> <flips|-> <out.json>   # build + integer lift + score
python3 g13_boolscan.py                          # 1156 single-boolean exact reduces (~22 min)
```
Modules: `gsym.py` (symbolic F_p forward eval), `gclose.py` (symbol closure),
`gred.py` (reduce_state: linear solve + substitution), `g14_print.py` (print the system).

## Next experiment
Post-process `boolscan.pkl` with the CORRECT criterion `ninc==0 and nzc==0 and residual==[]`
(the earlier "ZERO RESIDUAL at x2081" print was a false alarm — the obstruction moved into
the linear part, ninc=5). Then 2-bit and 3-bit combinations over the bits that reduce
`ninc + nzc + |residual|`. After that: make the booleans themselves symbolic (multilinear,
b^2=b) so the boolean search becomes algebra rather than enumeration.
