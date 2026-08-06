import heal_harness as H
import json
p=H.p
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F=H.fails()
print("39022 fails:", len(F), sorted(F))
# nonzero atoms
atoms=[]; reprs=[]; aeqs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line)
        atoms.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr','')); aeqs.append(dd.get('eqs',[]))
def atomval(items):
    s=0
    for m,c in items:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
nz=[i for i in range(len(atoms)) if atomval(atoms[i])!=0]
print("nonzero atoms at 39022:", len(nz), nz)
for i in nz:
    vs=set()
    for m,c in atoms[i]: vs.update(m)
    print(f"--- atom {i}: {reprs[i]}")
    print(f"     eqs={aeqs[i]}  vars={sorted(vs)}")
    for v in sorted(vs):
        print(f"       x_{v}={H.val[v]%p}  free={v in H.freeinp}  anc={sorted(H.anc.get(v,set()))[:6]}")
