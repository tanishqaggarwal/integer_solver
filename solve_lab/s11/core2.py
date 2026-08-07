import sys, os, json, itertools
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
def rep(v,tag):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]; F=L.failing_eqs(AV)
    print(f"{tag}: broken={B} failing={len(F)} score={L.NEQ-len(F)}")
    return B,F
base=load_raw(os.path.join(HERE,'data','modp5_out.json'))
rep(base,'modp5')
# the knobs feed the trio through gates; find, for each, the value that zeroes the target
def try_set(v, target_atom, knob):
    """binary/affine solve: response is affine in knob, so fit and invert"""
    f0=L.evalpoly(L.polys[target_atom],v)
    v1=list(v); L.ripple(v1,{knob:v[knob]+1}); f1=L.evalpoly(L.polys[target_atom],v1)
    d=f1-f0
    if d==0 or f0% d: return None
    return v[knob] - f0//d
best=None
for k688 in (7497,22820):
    for k1618 in (11436,14393):
        for k21617 in (14623,):
            v=list(base)
            s={}
            x=try_set(v,688,k688)
            if x is not None: L.ripple(v,{k688:x}); s['688']=k688
            x=try_set(v,1618,k1618)
            if x is not None: L.ripple(v,{k1618:x}); s['1618']=k1618
            x=lin_solve(21617,k21617,v)
            if x is not None: L.ripple(v,{k21617:x}); s['21617']=k21617
            B,F=rep(v,f"  knobs {k688}/{k1618}/{k21617} {s}")
            if best is None or len(F)<best[0]: best=(len(F),v,B)
n,v,B=best
print(f"BEST failing={n} score={L.NEQ-n} broken={B}")
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','core2_out.json'),'w'))
