"""Closure v3: repair a check by solving for ANY variable of the atom and then REALIZING that
value through the definer DAG down to a free input (depth-limited)."""
import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
CHK=[a for a in range(L.NA) if a not in L.atom_out]
FREE=set(u for u in range(L.NVARS) if u not in outs)
TP={u:i for i,u in enumerate(topo)}
ATP={a:max((TP[u] for u in L.avars[a]),default=0) for a in CHK}
CHK_ORD=sorted(CHK,key=lambda a:ATP[a])
def evalat(a,v):
    s=0
    for m,c in L.polys[a].items():
        t=c
        for w in m: t*=v[w]
        s+=t
    return s
def lin_in(a,u,v):
    c=0; rest=0
    for m,cc in L.polys[a].items():
        k=m.count(u)
        if k>1: return None
        t=cc
        if k==0:
            for w in m: t*=v[w]
            rest+=t
        else:
            for w in m:
                if w!=u: t*=v[w]
            c+=t
    return c,rest
def solve_for(a,u,v,target=0):
    r=lin_in(a,u,v)
    if r is None: return None
    c,rest=r
    if c==0: return None
    num=target-rest
    if num%c: return None
    return num//c
def realize(t,target,v,forbid,depth=3,trial=None):
    """Return dict of free-var assignments making v[t]==target, or None."""
    if trial is None: trial={}
    if t in FREE:
        if t in forbid: return None
        return {t:target}
    if depth<=0: return None
    d=outs.get(t)
    if d is None: return None
    old=v[t]; v[t]=target
    best=None
    for s in sorted(L.avars[d],key=lambda z:-TP[z]):
        if s==t: continue
        nv=solve_for(d,s,v,0)
        if nv is None or nv==v[s]: continue
        r=realize(s,nv,v,forbid,depth-1)
        if r is not None: best=r; break
    v[t]=old
    return best
def closure3(control, rounds=25, verbose=False, depth=3, forbid=()):
    ctrl=dict(control); forbid=set(forbid)|set(ctrl)
    v=[0]*L.NVARS
    for k,val in ctrl.items(): v[k]=val
    forward(v)
    best=None
    for r in range(rounds):
        changed=0
        for a in CHK_ORD:
            if evalat(a,v)==0: continue
            done=False
            # level 1: direct free variable
            cands=[]
            for u in L.avars[a]:
                if u not in FREE or u in forbid: continue
                nv=solve_for(a,u,v)
                if nv is not None and nv!=v[u]: cands.append((TP[u],u,nv))
            if cands:
                cands.sort(); _,u,nv=cands[-1]
                v[u]=nv; forward(v); changed+=1; done=True
            if done: continue
            # level 2+: solve for a gate output and realize it through its definer
            for u in sorted(L.avars[a],key=lambda z:-TP[z]):
                if u in FREE: continue
                nv=solve_for(a,u,v)
                if nv is None or nv==v[u]: continue
                r2=realize(u,nv,v,forbid,depth)
                if r2:
                    for k2,val2 in r2.items(): v[k2]=val2
                    forward(v)
                    if evalat(a,v)==0: changed+=1; break
        av=L.all_atom_values(v)
        nz=[a for a in CHK if av[a]!=0]
        sc=L.NEQ-len(L.failing_eqs(av))
        if verbose: print('  r%d score=%d nz=%d changed=%d'%(r,sc,len(nz),changed),flush=True)
        if best is None or sc>best[0]: best=(sc,list(v),list(nz))
        if not nz: return sc,v,nz
        if changed==0: break
    return best[0],best[1],best[2]
