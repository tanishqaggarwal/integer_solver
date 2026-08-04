import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
# apply core-reduction linear fix
H.val[14853]-=(H.val[14853]-H.val[12186])%p
H.val[16742]+=(H.val[24908]-H.val[16742])%p
H.forward()
V=H.val
F=set(H.fails())
print(f"after core fix: {len(F)} fails")
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
nz=[(i,atomval(i)) for i in range(len(atoms)) if atomval(i)!=0]
print(f"nonzero atoms: {len(nz)}")
for i,val in nz:
    a=atoms[i]
    deg=max(len(vs) for vs,_ in a['poly'])
    print(f"  atom{i} (deg{deg}, %p={'0' if val%p==0 else 'nz'}, n_eq={a['n_eq']}): {a['repr'][:75]}")
