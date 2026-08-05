import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
V=H.val
V[2081]=0
# also need to handle: x_6418,x_12553 now free (unloaded). zero the gadget-off handles.
H.forward()
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
def atomval(i):
    s=0
    for vs,c in atoms[i]['poly']:
        m=c
        for vi in vs: m*=V[vi]
        s+=m
    return s
nz=[i for i in range(len(atoms)) if atomval(i)!=0]
print(f"x_2081=0: nonzero atoms: {len(nz)}")
for i in nz:
    a=atoms[i]; val=atomval(i)
    print(f"  atom{i} (%p={'0' if val%p==0 else 'nz'}, n_eq={a['n_eq']}): {a['repr'][:85]}")
