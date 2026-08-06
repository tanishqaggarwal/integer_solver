#!/usr/bin/env python3
"""Drive the cascade directly: at each round, find nonzero atoms in failing eqs, and for each,
zero it by setting a FREE variable that appears with unit coefficient. Track the failing set —
does it plateau (bounded -> solvable) or truly explode?"""
import heal_harness as H, json
p=H.p
# load atoms
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
# map eq -> atoms
eq2atoms={}
for idx,d in enumerate(atoms):
    for e in d.get('eqs',[]): eq2atoms.setdefault(e,[]).append(idx)

def atomval(poly):
    s=0
    for vl,co in poly:
        t=co%p
        for x in vl: t=(t*H.val[x])%p
        s=(s+t)%p
    return s

vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[17325]=0; H.val[9413]=0; H.forward()
H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]
H.val[2964]=H.val[2099]; H.val[24548]=H.val[19964]
for s in (19569,11052):
    if s in H.freeinp: H.val[s]=0
H.forward()

everseen=set()
for rnd in range(40):
    F=set(H.fails()); everseen|=F
    if not F:
        print('*** ZERO FAILS ***'); break
    # collect nonzero atoms
    nzatoms=set()
    for e in F:
        for ai in eq2atoms.get(e,[]):
            if atomval(atoms[ai]['poly'])!=0: nzatoms.add(ai)
    fixed=0
    for ai in nzatoms:
        poly=atoms[ai]['poly']
        val=atomval(poly)
        if val==0: continue
        # find a free var with unit (±1) coefficient appearing alone (degree-1 term)
        knob=None; sign=None
        for vl,co in poly:
            if len(vl)==1 and vl[0] in H.freeinp and co%p in (1,p-1):
                knob=vl[0]; sign=1 if co%p==1 else -1; break
        if knob is None: continue
        # set knob so atom becomes 0: current val includes sign*knob; new_knob = knob - sign*val
        H.val[knob]=(H.val[knob]-sign*val)%p
        # keep it as a small-ish integer rep (signed)
        if H.val[knob]>p//2: H.val[knob]-=p
        fixed+=1
    H.forward()
    nf=len(H.fails())
    print(f'round {rnd}: fails={len(F)} nz_atoms={len(nzatoms)} fixed={fixed} -> new_fails={nf} everseen={len(everseen)}')
print(f'total distinct eqs ever in failing set (plateau size): {len(everseen)}')
