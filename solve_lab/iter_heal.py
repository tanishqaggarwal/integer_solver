import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('g1g2_closed.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
# load all atoms for residual tracking
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
def nzatoms():
    return [i for i in range(len(atoms)) if atomval(i)!=0]
print(f"start fails={len(H.fails())}, nz_atoms={len(nzatoms())}: {nzatoms()}")
# iterative residue heal: for each nonzero p-slack atom of form A - x_wire*x_partner (wire=p),
# set the free leaf's residue to match. Track fail count.
# General: find nonzero atoms, for each, try to zero via a free var in it (residue or exact).
for it in range(12):
    nz=nzatoms()
    F=set(H.fails())
    print(f"iter {it}: fails={len(F)}, nz_atoms={len(nz)}")
    if not nz: print("ALL ATOMS ZERO"); break
    changed=False
    for ai in nz:
        a=atoms[ai]; val=atomval(ai)
        # find a free leaf in this atom to adjust (residue mod p if p-granular partner available)
        vs_all=set()
        for vv,_ in a['poly']:
            vs_all|=set(vv)
        frees=[v for v in vs_all if v in H.freeinp]
        # try: for each free, does adjusting it by -val/coef zero the atom exactly?
        for fv in frees:
            base=atomval(ai)
            V[fv]+=1; fwd_from([fv]); slope=atomval(ai)-base; V[fv]-=1; fwd_from([fv])
            if slope!=0 and base%slope==0:
                V[fv]-=base//slope; fwd_from([fv]); changed=True; break
    if not changed:
        # try residue heal: pick a p-granular slack, match residues
        break
F=set(H.fails())
print(f"final: fails={len(F)}, nz_atoms={len(nzatoms())}")
