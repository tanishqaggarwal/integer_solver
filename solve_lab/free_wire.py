import json
p=2**256-2**32-977
with open('atoms/poly_atoms.jsonl') as f:
    atoms=[json.loads(l) for l in f]
from collections import defaultdict
eq_atoms=defaultdict(list)
for i,a in enumerate(atoms):
    for e in a['eqs']: eq_atoms[e].append(i)
# find atoms containing x_13859 and x_15616 (the slack multipliers)
for tgt in [13859,15616]:
    print(f"\n===== atoms containing x_{tgt} =====")
    for i,a in enumerate(atoms):
        for vs,c in a['poly']:
            if tgt in vs:
                deg=max(len(v) for v,_ in a['poly'])
                print(f"  atom{i} (deg{deg},n_eq={a['n_eq']}): {a['repr'][:90]}")
                break
