import json
import heal_harness as H
p=H.p
d=H.loadd('g1g2_closed.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
# atom footprints
with open('atoms/poly_atoms.jsonl') as f:
    atoms=[json.loads(l) for l in f]
def natom(vid):  # how many atoms contain vid
    c=0
    for a in atoms:
        for vs,_ in a['poly']:
            if vid in vs: c+=1; break
    return c
print("=== atom 7450: x_2964 - x_26756 - x_579 ===")
for n in [2964,26756,579]:
    print(f"  x_{n}: {'FREE' if n in H.freeinp else 'gate'} def={gdef.get(n,('',()))[0][:50]} n_atoms={natom(n)} val%p={V[n]%p!=0}")
print("=== atom 7452: 9367949*(x_24548 - x_25442) - x_7927 ===")
for n in [24548,25442,7927]:
    print(f"  x_{n}: {'FREE' if n in H.freeinp else 'gate'} def={gdef.get(n,('',()))[0][:50]} n_atoms={natom(n)}")
# residual values
print(f"\natom7450 = {V[2964]-V[26756]-V[579]}")
print(f"atom7452 = {9367949*(V[24548]-V[25442])-V[7927]}")
# verifier atoms full repr
print(f"\natom44342: {atoms[44342]['repr']}")
print(f"\natom45677: {atoms[45677]['repr']}")
