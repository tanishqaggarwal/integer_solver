"""Simultaneous exact repair: build the integer affine map free-vars -> atom residuals, solve over Z."""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
MUXV={22162,30213}
def mux(s,n=4):
    v=E.forward(s)
    for _ in range(n):
        s[22162]=v[13682]; s[30213]=v[18956]-v[32237]; v=E.forward(s)
    return v
def resid(s):
    s=dict(s); v=mux(s); return E.badatoms(v), v

def cone_free(aid):
    return set(E.cone(aid)[1])

def build(base, S, log=sys.stdout):
    r0,_=resid(base)
    cols={}; nonlin=set()
    t0=time.time()
    for k,f in enumerate(S):
        s1=dict(base); s1[f]=base.get(f,0)+1
        r1,_=resid(s1)
        s2=dict(base); s2[f]=base.get(f,0)+2
        r2,_=resid(s2)
        d1={}; d2={}
        for a in set(r0)|set(r1): d1[a]=r1.get(a,0)-r0.get(a,0)
        for a in set(r0)|set(r2): d2[a]=r2.get(a,0)-r0.get(a,0)
        for a in set(d1)|set(d2):
            if d2.get(a,0)!=2*d1.get(a,0): nonlin.add((f,a))
        cols[f]={a:x for a,x in d1.items() if x}
        if k%25==0: print(f"  probe {k}/{len(S)} {time.time()-t0:.0f}s",file=log,flush=True)
    return r0, cols, nonlin
