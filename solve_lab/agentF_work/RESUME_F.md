# RESUME_F — agent F (p-adic / multi-modular lifting)

## Best verified score: 39,026 / 39,033 (baseline, not yet improved)
Path: `solve_lab/best/new_instance_partial_39026.json` (shared baseline; re-verified by me).
No agentF assignment >= 39,026 produced yet.

## Established (mine, independently)
- Baseline claim CONFIRMED by running checker.py: 39026/39033, fails [12231,12270,12350,14584,18673,22044,29125].
- Own parser built: `agentF_work/parse.py|parse2.py|parse3.py`; artifacts `eqs3.pkl` (per-eq (mult,[(coef,atomAST)])).
- Instance = circuit. 39,033 eqs, each = scalar * core S; S = left-nested spine of <=26 atoms.
  96,883 distinct atoms in only 31 syntactic shapes (definitions x-f(..), booleans x(x-1), pins x-0/x-1).

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentF_work
    python3 parse3.py          # rebuilds eqs3.pkl (~30s)

## Next experiment
Build the definition DAG over atoms; find which variables are multiply-defined / cyclic; that is where
the residual must live. Then test the mod-p split independently (which p? derive from constants, do not assume).
