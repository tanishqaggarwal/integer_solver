# Agent F log — p-adic / multi-modular lifting angle

## 2026-08-07 start
- Verified baseline: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`. CONFIRMED.
- Built independent parser (`parse.py`,`parse2.py`,`parse3.py`) from raw text; no reliance on prior lab code.
  - 39,033 equations. Each LHS is a scalar multiple / square / (c1+c2) multiple of a single core S.
  - S is a left-nested spine `A0 + c1*A1 + c2*A2 + ...`, up to 26 atoms, small coefs.
  - 96,883 distinct atoms; 31 atom SHAPES only. Dominant shapes:
    (X-(X*X)) 160719, (X-(X+X)) 56409, (X-X) 45987, X 43773, (X-(C*X)) 28396,
    (X-(X-X)) 28082, (X-0) 20417, ((X*X)-X) 20242, (X-(1-X)) 17207, (X*X) 16413,
    (X*(1-X)) 14964, (X*(X-1)) 14950, (X-1) 13369, ...
  - So the instance is a straight-line circuit: most atoms are *definitions* x_out - f(inputs);
    plus boolean atoms X*(X-1) and constant pins (X-0),(X-1).
