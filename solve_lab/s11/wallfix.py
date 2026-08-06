import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
def lin_solve(a,t,v):
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc
        for u in m:
            if u!=t: term*=v[u]
        c+=term
    if c==0: return None
    old=v[t]; v[t]=0; rest=L.evalpoly(L.polys[a],v); v[t]=old
    if rest%c: return None
    return -rest//c
v=load_raw(os.path.join(HERE,'data','fix7_29539_7930.json'))
def rep(v,tag):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]; F=L.failing_eqs(AV)
    print(f"{tag}: broken={B} failing={len(F)} score={L.NEQ-len(F)}")
    for a in B: print(f"     a{a} in {len(L.atom2eq.get(a,{}))} eqs")
    return B,F
rep(v,'start')
for a,t in [(19297,30317),(19299,5146),(30984,2936),(21617,14623)]:
    x=lin_solve(a,t,v)
    print(f"  a{a} <- x{t}: {'SOLVABLE' if x is not None else 'NOT divisible'}")
    if x is not None: L.ripple(v,{t:x})
B,F=rep(v,'after wall handles')
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','wallfix_out.json'),'w'))
