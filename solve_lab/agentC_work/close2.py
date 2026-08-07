"""Closure v2: control seeds are fixed; all other free vars are re-solvable every round."""
import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
CHK=[a for a in range(L.NA) if a not in L.atom_out]
FREE=set(u for u in range(L.NVARS) if u not in outs)
TP={u:i for i,u in enumerate(topo)}
# topological rank of a check = max topo index of its inputs
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
def closure2(control, rounds=40, verbose=False, forbid=()):
    ctrl=dict(control); forbid=set(forbid)|set(ctrl)
    v=[0]*L.NVARS
    for k,val in ctrl.items(): v[k]=val
    forward(v)
    best=None
    for r in range(rounds):
        changed=0
        for a in CHK_ORD:
            val=evalat(a,v)
            if val==0: continue
            cands=[]
            for u in L.avars[a]:
                if u not in FREE or u in forbid: continue
                res=lin_in(a,u,v)
                if res is None: continue
                c,rest=res
                if c==0 or rest%c: continue
                nv=-rest//c
                if nv!=v[u]: cands.append((TP[u],u,nv))
            if not cands: continue
            cands.sort()
            _,u,nv=cands[-1]
            v[u]=nv; forward(v); changed+=1
        av=L.all_atom_values(v)
        nz=[a for a in CHK if av[a]!=0]
        sc=L.NEQ-len(L.failing_eqs(av))
        if verbose: print('  r%d score=%d nz=%d changed=%d'%(r,sc,len(nz),changed))
        if best is None or sc>best[0]: best=(sc,list(v),list(nz))
        if not nz: return sc,v,nz
        if changed==0: break
    return best[0],best[1],best[2]
