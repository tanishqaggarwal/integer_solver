# RESUME_J (agent J) — reduced-parameterization attack

## Status
- Verified baseline: `solve_lab/best/new_instance_partial_39026.json` = **39026/39033**
  (re-checked myself with solve_lab/checker.py; failing = [12231,12270,12350,14584,18673,22044,29125]).
- Best of my OWN: none better yet.
- Phase: independent re-parse of EQUATIONS.txt into atoms (jparse.py -> jmodel.pkl).

## Verdict on the claimed reduction ("thirteen 296-bit numbers")
- NOT YET EVALUATED.

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentJ_work
    python3 jparse.py          # ~60s, writes jmodel.pkl

## Next experiment
- Build definer graph from atoms, find free inputs, forward-evaluate; then test
  the claimed 13-parameter reduction by recomputing it.
