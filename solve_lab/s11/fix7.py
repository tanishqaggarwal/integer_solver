import sys, os, json, itertools, time
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
base=load_raw(os.path.join(HERE,'data','fix2_round.json'))
def report(v,tag):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]
    F=L.failing_eqs(AV)
    print(f"{tag}: broken atoms={len(B)} {B}  failing eqs={len(F)} score={L.NEQ-len(F)}")
    for a in B: print(f"      a{a} in {len(L.atom2eq.get(a,{}))} eqs")
    return len(F)
report(base,'fix2_round')
for combo in [(29539,),(7930,),(29539,7930)]:
    v=list(base)
    seeds={}
    for a in combo:
        t = 14853 if a==29539 else 24548
        x=lin_solve(a,t,v)
        seeds[t]=x
    L.ripple(v,seeds)
    n=report(v,f"  fix {combo} via free mirrors")
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)},
              open(os.path.join(HERE,'data','fix7_%s.json'%'_'.join(map(str,combo))),'w'))
