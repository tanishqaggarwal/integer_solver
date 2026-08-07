"""Closure v4: supports DETACHING gate variables (their definer becomes a broken check and the
variable becomes a controllable free input)."""
import sys, json, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
TP={u:i for i,u in enumerate(topo)}
ALL=list(range(L.NA))
def make_frame(detach):
    det=set(detach)
    outs2={x:a for x,a in outs.items() if x not in det}
    lin2={x:LN[x] for x in LN if x not in det}
    gorder=[u for u in topo if u in lin2]
    chk=[a for a in ALL if a not in L.atom_out or L.atom_out[a][1] in det]
    chk.sort(key=lambda a: max((TP[u] for u in L.avars[a]),default=0))
    return outs2,lin2,gorder,chk
def mk(detach=()):
    outs2,lin2,gorder,chk=make_frame(detach)
    FREE=set(u for u in range(L.NVARS) if u not in outs2)
    def fwd(v):
        for u in gorder:
            c,rest=lin2[u]; s=0
            for cc,m in rest:
                t=cc
                for w in m: t*=v[w]
                s+=t
            v[u]=-s//c if s%c==0 else 0
        return v
    return outs2,FREE,fwd,chk
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
_v0=[0]*L.NVARS
_o2,_F,_fw,_C=mk(())
_fw(_v0)
WIRES={u for u in range(L.NVARS) if _v0[u]==P}
def closure4(control, detach=(), rounds=25, verbose=False, depth=3, extra_forbid=()):
    outs2,FREE,fwd,CHK=mk(detach)
    forbid=set(control)|set(detach)|WIRES|set(extra_forbid)
    def realize(t,target,v,depth):
        if t in FREE:
            return None if t in forbid else {t:target}
        if depth<=0: return None
        d=outs2.get(t)
        if d is None: return None
        old=v[t]; v[t]=target; best=None
        for s in sorted(L.avars[d],key=lambda z:-TP[z]):
            if s==t or s in WIRES: continue
            nv=solve_for(d,s,v,0)
            if nv is None or nv==v[s]: continue
            r=realize(s,nv,v,depth-1)
            if r is not None: best=r; break
        v[t]=old
        return best
    v=[0]*L.NVARS
    for k,val in control.items(): v[k]=val
    fwd(v)
    best=None
    for r in range(rounds):
        changed=0
        for a in CHK:
            if evalat(a,v)==0: continue
            cands=[]
            for u in L.avars[a]:
                if u not in FREE or u in forbid: continue
                nv=solve_for(a,u,v)
                if nv is not None and nv!=v[u]: cands.append((TP[u],u,nv))
            if cands:
                cands.sort(); _,u,nv=cands[-1]
                v[u]=nv; fwd(v); changed+=1; continue
            for u in sorted(L.avars[a],key=lambda z:-TP[z]):
                if u in FREE: continue
                nv=solve_for(a,u,v)
                if nv is None or nv==v[u]: continue
                r2=realize(u,nv,v,depth)
                if r2:
                    for k2,val2 in r2.items(): v[k2]=val2
                    fwd(v)
                    if evalat(a,v)==0: changed+=1; break
        av=L.all_atom_values(v)
        nz=[a for a in range(L.NA) if av[a]!=0]
        sc=L.NEQ-len(L.failing_eqs(av))
        if verbose: print('  r%d score=%d nz=%d changed=%d'%(r,sc,len(nz),changed),flush=True)
        if best is None or sc>best[0]: best=(sc,list(v),list(nz))
        if not nz: return sc,v,nz
        if changed==0: break
    return best[0],best[1],best[2]
