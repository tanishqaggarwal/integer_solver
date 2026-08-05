import heal_harness as H
import json
p=H.p
ATOMS=[]; reprs=[]; aeqs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr','')); aeqs.append(dd.get('eqs',[]))
d=H.loadd('best/new_instance_partial_39021.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
print("39021 fails:", len(H.fails()), sorted(H.fails()))
for a in [30278,30280,42405]:
    print(f"\n=== atom {a}: {reprs[a]}")
    print(f"   eqs={aeqs[a]}")
    vs=sorted(set(x for m,c in ATOMS[a] for x in m))
    for v in vs:
        print(f"   x_{v}={H.val[v]%p}  free={v in H.freeinp}  anc={sorted(H.anc.get(v,set()))[:6]}")
# check x_24601 and the pin targets
print("\nx_24601 (control bit) =", H.val[24601])
