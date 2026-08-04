import heal_harness as H
import json, random
from collections import defaultdict
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()

# load atoms
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line)
        poly=[(tuple(m),c) for m,c in dd['poly']]
        atoms.append(poly)
print("natoms:", len(atoms))

def atomval(poly,val):
    s=0
    for m,c in poly:
        t=c%p
        for v in m: t=(t*val[v])%p
        s=(s+t)%p
    return s

# nonzero atoms at 39013
nz=[i for i,a in enumerate(atoms) if atomval(a,H.val)!=0]
print("nonzero atoms at 39013:", len(nz), nz[:50])

# breakable atoms via random free perturbation
breakable=set()
freelist=sorted(H.freeinp)
for seed in range(4):
    random.seed(100+seed)
    for v in freelist: H.val[v]=random.randrange(p)
    H.forward()
    for i,a in enumerate(atoms):
        if atomval(a,H.val)!=0: breakable.add(i)
# restore
for v in H.freeinp: H.val[v]=base[v]
H.forward()
print("breakable atoms:", len(breakable))
print("nonzero subset breakable:", set(nz)<=breakable)
import pickle
pickle.dump({'atoms':atoms,'nz':nz,'breakable':sorted(breakable)},
    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atoms.pkl','wb'))
