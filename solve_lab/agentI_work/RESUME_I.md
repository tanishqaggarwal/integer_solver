# Agent I — RESUME (build-from-scratch / complete-search angle)

## Status
- Baseline re-verified myself: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> 39026/39033, failing [12231, 12270, 12350, 14584, 18673, 22044, 29125]. CONFIRMED.
- Best of my own so far: none above baseline yet.

## My independent model (built from scratch, do not redo)
- `agentI_work/parse.py` -> `atoms.pkl`: 39,033 eqs, **40,885 distinct atoms**.
  Each equation = outer wrapper (c*, P*P, c1*P+c2*P) applied to a CORE = left-nested
  sum `a0 + c1*a1 + c2*a2 + ...` of atoms.
- `agentI_work/poly.py` -> `polys.pkl`: every atom is a polynomial of **degree <= 2 with
  <= 3 terms**. Histogram (deg,nterms): (1,1)2616 (1,2)8952 (1,3)9067 (2,1)735
  (2,2)19003 (2,3)512.
- Equation atom-count: min 1, mostly 3-24.

## Commands to re-enter
```
cd /home/user/integer_solver/solve_lab/agentI_work
python3 parse.py      # ~35s, writes atoms.pkl
python3 poly.py       # writes polys.pkl
python3 dag.py        # gate DAG / free inputs
```

## Next experiment
Build the gate DAG, identify free inputs, run full propagation from empty, and locate
the exact constraint core that blocks. Then complete search (CDCL/no-good) over that core.
