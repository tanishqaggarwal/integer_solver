"""Batch fixed-point repair.

Every surviving broken atom is a check of the form 'this FREE variable must equal that computed
value'.  Repairing them one at a time cascades forever because the copies form a network.  So do
it Jacobi-style: in each pass compute, for every broken check, the value its free variable must
take, apply them ALL as seeds, ripple the gates once, and repeat.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)

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
    old=v[t]; v[t]=0
    rest=L.evalpoly(L.polys[a],v)
    v[t]=old
    if rest%c: return None
    return -rest//c

def broken(v):
    return [a for a in range(L.NA) if L.evalpoly(L.polys[a],v)!=0]

def passes(v, npass=40, verbose=True):
    hist=[]
    for it in range(npass):
        B=broken(v)
        if verbose: print(f"  pass{it}: broken={len(B)} {B[:14]}", flush=True)
        if not B: return v,[]
        if len(hist)>=6 and len(set(hist[-6:]))==1 and hist[-1]==len(B):
            print("  stalled"); return v,B
        hist.append(len(B))
        seeds={}
        for a in B:
            cands=[t for t in sorted(L.avars[a]) if t in FREE and t not in seeds]
            got=False
            for t in cands:
                x=lin_solve(a,t,v)
                if x is not None and x!=v[t]:
                    seeds[t]=x; got=True; break
            if not got:
                # fall back to the variable this atom defines
                o=L.atom_out.get(a)
                if o is not None:
                    t=o[1]
                    x=lin_solve(a,t,v)
                    if x is not None and x!=v[t] and t not in seeds: seeds[t]=x
        if not seeds:
            print("  no seeds"); return v,B
        L.ripple(v, seeds)
    return v, broken(v)

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'data','fix2_round.json')
    v=load_raw(src)
    t0=time.time()
    v,B=passes(v)
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    F=L.failing_eqs(AV)
    print(f"RESULT broken atoms={len(B)} failing eqs={len(F)} score={L.NEQ-len(F)} ({time.time()-t0:.0f}s)")
    out=os.path.join(HERE,'data','fixpoint_out.json')
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(out,'w')); print("saved",out)
