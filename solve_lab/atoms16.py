import json
from collections import defaultdict
p=2**256-2**32-977
d=json.load(open('g1g2_closed.json'))
val=defaultdict(int)
for k,vv in d.items(): val[int(k[2:])]=int(vv)
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
def ev(poly):
    s=0
    for vs,c in poly:
        m=c
        for vi in vs: m*=val[vi]
        s+=m
    return s
nz=[]
for i,a in enumerate(atoms):
    v=ev(a['poly'])
    if v!=0: nz.append((i,v,a))
print(f"nonzero atoms: {len(nz)}")
for i,v,a in nz:
    print(f"\natom {i}: {a['repr'][:110]}")
    print(f"   %p={'0' if v%p==0 else 'nonzero'}, n_eq={a['n_eq']}, eqs={a['eqs'][:8]}")
