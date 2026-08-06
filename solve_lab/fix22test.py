import heal_harness as H
import json
p=H.p
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F0=set(H.fails())
print("before:", len(F0), "fails")
tgt2099=H.val[2099]; tgt19964=H.val[19964]
# set free inputs x_7068 := x_2099, x_4432 := x_19964
H.val[7068]=tgt2099
H.val[4432]=tgt19964
H.forward()
F1=set(H.fails())
print("after G1/G2 fix:", len(F1), "fails")
print("newly broken:", sorted(F1-F0))
print("newly fixed:", sorted(F0-F1))
# nonzero atoms now
atoms=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); atoms.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
def atomval(items):
    s=0
    for m,c in items:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
nz=[i for i in range(len(atoms)) if atomval(atoms[i])!=0]
print("nonzero atoms after fix:", len(nz), nz)
# are the new-broken atoms mod-p or integer-carry?
for i in nz:
    print(f"  atom {i}: {reprs[i][:80]}")
