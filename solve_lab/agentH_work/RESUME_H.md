# Agent H — RESUME (decomposition angle)

## Best verified score
39,026 / 39,033 — the inherited deliverable, re-verified by me this session:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> satisfied 39026/39033, failing [12231, 12270, 12350, 14584, 18673, 22044, 29125].
I have not yet beaten it. No agentH artifact yet exceeds it.

## Established (my own computation, independent of prior sessions)
- Parsed EQUATIONS.txt -> 42,267 atoms / 39,033 equations (`agentH_work/model.py`, cache model.pkl).
  Matches prior sessions' atom count.
- Eq-var bipartite graph: ONE connected component over all 38,748 vars. Naive component
  decomposition is worthless. (`graph1.py`)
- Atom-eq graph: 1 giant component (39,033 eqs) + 3,234 singleton atoms (the single-atom
  equations that force an atom to 0).
- Atom degrees: 19,780 linear, 21,788 quadratic, 699 quartic. 10,478 equations are squares.
- 36,096 atoms have >=1 "definition candidate" (a var with coeff +-1 not occurring elsewhere in
  the atom); 6,171 atoms have none (pure checks).

## Re-entry
cd /home/user/integer_solver/solve_lab/agentH_work && python3 model.py   # rebuilds caches
Caches: model.pkl (atoms+eqs), polys.pkl (atom polynomials), defcands.pkl.

## Next experiment
Build the gate DAG (choose one definition per defined var), identify free inputs, then compute
the free-input hypergraph induced by the CHECK atoms and look for small separators / low
treewidth blocks around the 7-equation residual.
