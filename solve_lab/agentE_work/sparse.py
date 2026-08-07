"""Sparse exact integer solve: singleton + unit-pivot elimination with column index, HNF core."""
import sys, time, heapq
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import intsolve

def solve_sparse(rows, rhs, names=None, verbose=True, maxcore=400, maxbits=200000):
    R=[dict(r) for r in rows]; B=list(rhs)
    n=len(R); alive=set(i for i in range(n))
    colidx={}
    for i,r in enumerate(R):
        for v in r: colidx.setdefault(v,set()).add(i)
    fixed={}; subs=[]
    heap=[(len(R[i]),i) for i in alive]; heapq.heapify(heap)
    t0=time.time(); nsing=0
    def kill(i):
        alive.discard(i)
        for v in R[i]: colidx[v].discard(i)
    def setvar(v,val):
        for j in list(colidx.get(v,())):
            if j in alive and v in R[j]:
                c=R[j].pop(v); B[j]-=c*val
                heapq.heappush(heap,(len(R[j]),j))
        colidx[v]=set()
    while True:
        # process rows by ascending size
        prog=False
        while heap:
            sz,i=heapq.heappop(heap)
            if i not in alive: continue
            if len(R[i])!=sz:
                heapq.heappush(heap,(len(R[i]),i)); 
                if len(R[i])>sz: continue
            if not R[i]:
                if B[i]!=0: return None,'row %s unsatisfiable'%(names[i] if names else i),fixed
                kill(i); prog=True; continue
            if len(R[i])==1:
                v,c=next(iter(R[i].items()))
                if B[i]%c: return None,'row %s: rhs %% %s != 0'%(names[i] if names else i,c),fixed
                val=B[i]//c
                kill(i); fixed[v]=val; nsing+=1; prog=True
                setvar(v,val); continue
            # unit pivot in this row?
            pv=None
            for v,c in R[i].items():
                if c==1 or c==-1: pv=(v,c); break
            if pv is None:
                continue   # leave for core
            v,c=pv
            rest={u:x for u,x in R[i].items() if u!=v}; rr=B[i]
            subs.append((v,rest,rr,c)); kill(i)
            for j in list(colidx.get(v,())):
                if j not in alive or v not in R[j]: continue
                f=R[j].pop(v); q=f//c
                for u,x in rest.items():
                    nv=R[j].get(u,0)-q*x
                    if nv==0:
                        if u in R[j]: del R[j][u]; colidx[u].discard(j)
                    else:
                        if u not in R[j]: colidx.setdefault(u,set()).add(j)
                        R[j][u]=nv
                B[j]-=q*rr
                if B[j] and B[j].bit_length() > maxbits:
                    return None,'coefficient blowup',fixed
                heapq.heappush(heap,(len(R[j]),j))
            colidx[v]=set()
            prog=True
        if not prog: break
        heap=[(len(R[i]),i) for i in alive]; heapq.heapify(heap)
    core=[i for i in alive if R[i]]
    corevars=sorted(set().union(*[set(R[i]) for i in core])) if core else []
    if verbose: print(f"  elim: {nsing} singletons, {len(subs)} pivots, core {len(core)}x{len(corevars)}, {time.time()-t0:.1f}s",flush=True)
    if len(core)>maxcore or len(corevars)>maxcore:
        return None,'core too large (%dx%d)'%(len(core),len(corevars)),fixed
    if core:
        mb=max((abs(x).bit_length() for i in core for x in R[i].values()), default=0)
        mb=max(mb, max((abs(B[i]).bit_length() for i in core), default=0))
        if mb>4000:
            return None,'core entries too large (%d bits)'%mb,fixed
    sol=dict(fixed)
    if core:
        A=[[R[i].get(v,0) for v in corevars] for i in core]
        b=[B[i] for i in core]
        d,ker=intsolve.solve_int(A,b)
        if d is None: return None,'core infeasible',fixed
        for v,x in zip(corevars,d): sol[v]=int(x)
    for v,rest,rr,c in reversed(subs):
        s=rr-sum(x*sol.get(u,0) for u,x in rest.items())
        if s%c: return None,'backsub non-integral',sol
        sol[v]=s//c
    return sol,'ok',sol
