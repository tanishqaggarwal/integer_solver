"""Sparse exact integer solve: singleton + unit-pivot elimination, HNF on the residual core."""
import sys, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import intsolve

def solve_sparse(rows, rhs, names=None, verbose=True):
    R=[dict(r) for r in rows]; B=list(rhs)
    n=len(R); alive=set(range(n))
    fixed={}      # var -> forced value
    subs=[]       # (var, rest, rhs, coef)
    def substitute(v, val):
        for j in list(alive):
            if v in R[j]:
                c=R[j].pop(v); B[j]-=c*val
    t0=time.time(); nsing=0
    while True:
        prog=False
        for i in list(alive):
            if i not in alive: continue
            if not R[i]:
                if B[i]!=0:
                    return None, 'row %s unsatisfiable'%(names[i] if names else i), fixed
                alive.discard(i); prog=True; continue
            if len(R[i])==1:
                v,c=next(iter(R[i].items()))
                if B[i]%c:
                    return None, 'row %s: %s not divisible by %s'%(names[i] if names else i, B[i], c), fixed
                val=B[i]//c
                alive.discard(i); fixed[v]=val; nsing+=1; prog=True
                substitute(v,val)
        if prog: continue
        best=None
        for i in alive:
            for v,c in R[i].items():
                if abs(c)==1:
                    k=len(R[i])
                    if best is None or k<best[0]: best=(k,i,v,c)
                    break
        if best is None: break
        _,i,v,c=best; rest={u:x for u,x in R[i].items() if u!=v}; rr=B[i]
        subs.append((v,rest,rr,c)); alive.discard(i)
        for j in list(alive):
            if v in R[j]:
                f=R[j].pop(v); q=f//c
                for u,x in rest.items():
                    R[j][u]=R[j].get(u,0)-q*x
                    if R[j][u]==0: del R[j][u]
                B[j]-=q*rr
    core=[i for i in alive if R[i]]
    if verbose: print(f"  elim: {nsing} singletons, {len(subs)} unit pivots, core {len(core)} rows, {time.time()-t0:.1f}s",flush=True)
    corevars=sorted(set().union(*[set(R[i]) for i in core])) if core else []
    sol=dict(fixed)
    if core:
        if verbose: print(f"  core vars {len(corevars)}",flush=True)
        A=[[R[i].get(v,0) for v in corevars] for i in core]
        b=[B[i] for i in core]
        t1=time.time()
        d,ker=intsolve.solve_int(A,b)
        if verbose: print(f"  core HNF {time.time()-t1:.1f}s feasible={d is not None}",flush=True)
        if d is None: return None,'core infeasible',fixed
        for v,x in zip(corevars,d): sol[v]=int(x)
    for v,rest,rr,c in reversed(subs):
        s=rr-sum(x*sol.get(u,0) for u,x in rest.items())
        if s%c: return None,'backsub non-integral',sol
        sol[v]=s//c
    return sol,'ok',sol
