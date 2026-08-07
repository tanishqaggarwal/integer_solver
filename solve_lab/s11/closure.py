"""How big is the CLOSED region (equations closed under the variables they contain)?"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
AV=[atomval(a,v) for a in range(L.NA)]
FAIL=[e for e in range(L.NEQ) if sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())!=0]
R=set(FAIL); seenv=set()
t0=time.time()
for it in range(30):
    vs=set()
    for e in R:
        for a in L.eq_atoms[e][2]: vs |= set(L.avars[a])
    new=vs-seenv; seenv|=vs
    add=set()
    for u in new:
        for a in L.var_atoms[u]: add |= set(L.atom2eq.get(a,{}))
    before=len(R); R|=add
    print(f"  it{it}: eqs={len(R)} vars={len(seenv)} (+{len(R)-before} eqs) ({time.time()-t0:.0f}s)", flush=True)
    if len(R)==before: break
print(f"CLOSURE: {len(R)} of {L.NEQ} equations, {len(seenv)} of {L.NVARS} variables")
