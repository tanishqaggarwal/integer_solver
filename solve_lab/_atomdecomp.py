import heal_harness as H, json
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
# load atoms with their eqs membership
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        atoms.append(json.loads(line))
# for each failing eq, find atoms that list it in eqs
eq_atoms=defaultdict={}
from collections import defaultdict
eq_atoms=defaultdict(list)
for ai,a in enumerate(atoms):
    for e in a['eqs']:
        if e in F: eq_atoms[e].append(ai)
def evalpoly(poly):
    s=0
    for mon,coef in poly:
        t=coef
        for vv in mon: t*=H.val[vv]
        s+=t
    return s
for e in F:
    print(f"\n=== eq {e}: {len(eq_atoms[e])} atoms ===")
    nz=[]
    for ai in eq_atoms[e]:
        a=atoms[ai]; val=evalpoly(a['poly'])
        deg=max((len(m) for m,_ in a['poly']),default=0)
        if val%p!=0:
            nz.append((ai,deg,a['repr'][:60]))
    print(f"  nonzero atoms(mod p): {len(nz)}")
    for ai,deg,r in nz: print(f"    atom {ai} deg{deg}: {r}")
