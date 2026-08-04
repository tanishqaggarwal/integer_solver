import json
from collections import defaultdict
p=2**256-2**32-977
d=json.load(open('sy_regime11_39018.json'))
val=defaultdict(int)
for k,vv in d.items(): val[int(k[2:])]=int(vv)
# load atoms
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
def ev(poly):
    s=0
    for term in poly:
        vs,c=term
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
    print(f"\natom {i}: {a['repr'][:120]}")
    print(f"   value%p = {v%p}")
    print(f"   value: {'p-mult' if v%p==0 else 'sub-p'}, n_eq={a['n_eq']}, eqs={a['eqs'][:8]}")
