"""How much does the 'compensator disturbs nothing outside the region' filter cost?
   Count candidates rejected by it, and how far the region would have to grow."""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import atomval, load_raw, deltas
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
LAB=os.path.join(HERE,'..')
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
AV=[atomval(a,v) for a in range(L.NA)]
def eqs(e): return sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())
FAIL=[e for e in range(L.NEQ) if eqs(e)!=0]
core=set()
for e in FAIL:
    for a in L.eq_atoms[e][2]: core |= set(L.avars[a])
collat=set()
uc=[]
for u in sorted(core):
    d1=deltas(v,AV,u,1)
    if not d1: continue
    d2=deltas(v,AV,u,2)
    if not all(d2.get(e,0)==2*d1.get(e,0) for e in set(d1)|set(d2)): continue
    uc.append(u); collat |= {e for e in d1 if e not in set(FAIL)}
comp=set()
for e in collat:
    for a in L.eq_atoms[e][2]: comp |= set(L.avars[a])
comp -= set(uc)
kept=0; rejected=0; extra=collections.Counter()
for u in sorted(comp):
    d1=deltas(v,AV,u,1)
    if not d1 or not (set(d1)&collat): continue
    d2=deltas(v,AV,u,2)
    if not all(d2.get(e,0)==2*d1.get(e,0) for e in set(d1)|set(d2)): continue
    out=[e for e in d1 if e not in collat and e not in set(FAIL)]
    if out: rejected+=1; extra[len(out)]+=1
    else: kept+=1
print(f"core exact-linear vars : {len(uc)}   collateral eqs : {len(collat)}")
print(f"compensators kept (no new collateral) : {kept}")
print(f"compensators REJECTED by that filter  : {rejected}")
print(f"  extra-equation counts they'd add    : {dict(sorted(extra.items())[:10])}")
