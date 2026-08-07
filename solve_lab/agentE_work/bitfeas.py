"""Per-bit pin-system feasibility: is the linear Diophantine repair of one activated bit solvable?"""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast, sparse
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
CORE={20212,20215,24403,747,4872,4877}
def cone_free(a): return set(E.cone(a)[1])
def closure(base_seed, frozen, maxvars=6000):
    v0=E.forward(base_seed); bad0=E.badatoms(v0)
    S=[]; cols={}; nonlin=set(); processed=set()
    pending=set(bad0)-CORE
    for rnd in range(25):
        newS=set()
        for a in pending: newS|=cone_free(a)
        newS-=set(S)|frozen
        newS=sorted(newS)
        if not newS: break
        for f in newS:
            b1,_=fast.resid_delta(v0,bad0,{f:v0[f]+1})
            b2,_=fast.resid_delta(v0,bad0,{f:v0[f]+2})
            col={}
            for a in set(b1)|set(bad0):
                d=b1.get(a,0)-bad0.get(a,0)
                if d: col[a]=d
            for a in set(b2)|set(bad0)|set(col):
                if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0): nonlin.add((f,a))
            cols[f]=col; S.append(f)
        aff=set()
        for f in newS: aff|=set(cols[f])
        processed|=pending
        pending=(aff|set(bad0))-processed-CORE
        if len(S)>maxvars: break
    return v0,bad0,S,cols,nonlin
def feas(bit, log=sys.stdout):
    base={18956:C, bit:1}
    v0,bad0,S,cols,nonlin=closure(base,{18956,bit})
    Sset=set(S); atoms=set(bad0)-CORE
    for f in S: atoms|=set(cols[f])
    atoms-=CORE
    nl={a for f,a in nonlin if f in Sset}
    use=sorted(a for a in atoms if a not in nl)
    rowmap={a:{} for a in use}
    for f in S:
        for a,c in cols[f].items():
            if a in rowmap: rowmap[a][f]=c
    rows=[rowmap[a] for a in use]; rhs=[-bad0.get(a,0) for a in use]
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=False)
    return sol,msg,S,use,bad0
if __name__=='__main__':
    bits=[int(x) for x in sys.argv[1:]]
    for b in bits:
        t0=time.time()
        sol,msg,S,use,bad0=feas(b)
        nz=sum(1 for v in (sol or {}).values() if v)
        print(f"bit x_{b}: vars={len(S)} rows={len(use)} -> {msg[:110]} moved={nz} ({time.time()-t0:.0f}s)",flush=True)
