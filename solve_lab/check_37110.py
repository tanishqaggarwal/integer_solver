import json
p=2**256-2**32-977
with open('atoms/poly_atoms.jsonl') as f:
    atoms=[json.loads(l) for l in f]
# equations containing atom 37110
eqs37110=atoms[37110]['eqs']
print(f"atom 37110 (x_26064 - p) in eqs: {eqs37110}")
# for each such equation, find its atom decomposition. Need eq->atoms map.
# build from atoms' eqs lists
eq_atoms=defaultdict={} 
from collections import defaultdict
eq_atoms=defaultdict(list)
for i,a in enumerate(atoms):
    for e in a['eqs']: eq_atoms[e].append(i)
# examine eq 8429
for E in eqs37110[:3]:
    ats=eq_atoms[E]
    print(f"\n=== eq {E}: {len(ats)} atoms ===")
    for ai in ats:
        a=atoms[ai]
        # is it a product-slack (has degree-2 term) with a potentially-free partner?
        deg=max(len(vs) for vs,_ in a['poly'])
        print(f"  atom{ai} (deg{deg}): {a['repr'][:80]}")
