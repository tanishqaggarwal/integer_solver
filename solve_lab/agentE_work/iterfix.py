"""Iterated closure+solve (Newton over the exact linear model) from an arbitrary seed."""
import sys, time, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B, sparse
C=B.C
def iterate(base, frozen, iters=6, maxr=5, maxv=4000, exclude=None, log=sys.stdout):
    exclude=exclude or set()
    s=dict(base); hist=[]
    for it in range(iters):
        v0=E.forward(s); bad0=E.badatoms(v0)
        bad0={a:r for a,r in bad0.items() if a not in exclude}
        ff=E.eqfails(E.badatoms(v0))
        print(f"  it{it}: bad={sorted(bad0)} fails={len(ff)}",file=log,flush=True)
        hist.append((sorted(bad0),len(ff)))
        if not bad0: return s,hist,True
        S,cols,nonlin,rounds=B.build_closure_from(v0,bad0,frozen,maxr,maxv) if hasattr(B,'build_closure_from') else (None,)*4
        if S is None:
            # inline closure
            import fast
            S=[]; cols={}; nonlin=set(); processed=set(); rounds={}
            pending=set(bad0)
            for rnd in range(maxr+1):
                newS=set()
                for a in pending: newS|=set(E.cone(a)[1])
                newS-=set(S)|frozen
                newS=sorted(newS)
                if not newS: break
                for f in newS:
                    b1,_=fast.resid_delta(v0,bad0,{f:v0[f]+1})
                    b2,_=fast.resid_delta(v0,bad0,{f:v0[f]+2})
                    col={}
                    for a in set(b1)|set(bad0):
                        d=b1.get(a,0)-bad0.get(a,0)
                        if d and a not in exclude: col[a]=d
                    for a in set(b2)|set(bad0)|set(col):
                        if a in exclude: continue
                        if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0): nonlin.add((f,a))
                    cols[f]=col; S.append(f); rounds.setdefault(rnd,[]).append(f)
                aff=set()
                for f in newS: aff|=set(cols[f])
                processed|=pending; pending=(aff|set(bad0))-processed
                if len(S)>maxv: break
        got=None
        for mr in sorted(rounds):
            Sp=[]
            for r in sorted(rounds):
                if r<=mr: Sp+=rounds[r]
            Sset=set(Sp); atoms=set(bad0)
            for f in Sp: atoms|=set(cols[f])
            nl={a for f,a in nonlin if f in Sset}
            use=sorted(a for a in atoms if a not in nl and a not in exclude)
            rowmap={a:{} for a in use}
            for f in Sp:
                for a,c in cols[f].items():
                    if a in rowmap: rowmap[a][f]=c
            rows=[rowmap[a] for a in use]; rhs=[-bad0.get(a,0) for a in use]
            sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=False)
            if sol is None: continue
            ns=dict(s)
            for f,d in sol.items():
                if d: ns[f]=v0[f]+d
            vv=E.forward(ns); aav=E.badatoms(vv)
            aav={a:r for a,r in aav.items() if a not in exclude}
            nf=len(E.eqfails(E.badatoms(vv)))
            if got is None or (len(aav),nf)<(len(got[2]),got[3]): got=(msg,ns,aav,nf,mr)
            if not aav: break
        if got is None:
            print(f"  it{it}: no linear solution",file=log,flush=True)
            return s,hist,False
        print(f"  it{it}: r<={got[4]} -> bad={sorted(got[2])} fails={got[3]}",file=log,flush=True)
        if sorted(got[2])==sorted(bad0):
            return got[1],hist,False
        s=got[1]
    return s,hist,False
