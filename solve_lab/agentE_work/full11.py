"""Joint (1,1) solve: activate one a-tree bit and one b-tree bit, close the MUX, then solve the
   whole pin system as one exact linear Diophantine system over Z.  No atoms excluded."""
import sys, time, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast, sparse
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
NOCORE=set()

def mux11(s,n=6):
    v=E.forward(s)
    for _ in range(n):
        s[22162]=v[13682]; s[30213]=v[18956]-v[32237]; v=E.forward(s)
    return v,s

def closure(v0,bad0,frozen,maxr=6,maxv=6000,verbose=False):
    S=[]; cols={}; nonlin=set(); processed=set(); rounds={}
    pending=set(bad0)
    for rnd in range(maxr+1):
        newS=set()
        for a in pending: newS|=set(E.cone(a)[1])
        newS-=set(S)|frozen
        newS=sorted(newS)
        if not newS: break
        t0=time.time()
        for f in newS:
            b1,_=fast.resid_delta(v0,bad0,{f:v0[f]+1})
            b2,_=fast.resid_delta(v0,bad0,{f:v0[f]+2})
            col={}
            for a in set(b1)|set(bad0):
                d=b1.get(a,0)-bad0.get(a,0)
                if d: col[a]=d
            for a in set(b2)|set(bad0)|set(col):
                if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0): nonlin.add((f,a))
            cols[f]=col; S.append(f); rounds.setdefault(rnd,[]).append(f)
        aff=set()
        for f in newS: aff|=set(cols[f])
        processed|=pending; pending=(aff|set(bad0))-processed
        if verbose: print(f"   rnd{rnd}: +{len(newS)} ({len(S)}) {time.time()-t0:.0f}s",flush=True)
        if len(S)>maxv: break
    return S,cols,nonlin,rounds

def solve_pair(abit,bbit,verbose=True):
    s={18956:C, abit:1, bbit:1}
    v0,s=mux11(s)
    bad0=E.badatoms(v0)
    if verbose: print(f"[{abit},{bbit}] start bad={sorted(bad0)} fails={len(E.eqfails(bad0))}",flush=True)
    S,cols,nonlin,rounds=closure(v0,bad0,{18956,abit,bbit},verbose=verbose)
    best=None
    for maxr in sorted(rounds):
        Sp=[]
        for r in sorted(rounds):
            if r<=maxr: Sp+=rounds[r]
        Sset=set(Sp); atoms=set(bad0)
        for f in Sp: atoms|=set(cols[f])
        nl={a for f,a in nonlin if f in Sset}
        use=sorted(a for a in atoms if a not in nl)
        rowmap={a:{} for a in use}
        for f in Sp:
            for a,c in cols[f].items():
                if a in rowmap: rowmap[a][f]=c
        rows=[rowmap[a] for a in use]; rhs=[-bad0.get(a,0) for a in use]
        sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=False)
        if verbose: print(f"   r<={maxr}: vars={len(Sp)} rows={len(use)} nl={len(nl)} -> {msg[:80]}",flush=True)
        if sol is None: continue
        ns=dict(s)
        for f,d in sol.items():
            if d: ns[f]=v0[f]+d
        v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
        if verbose: print(f"   EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:12]}",flush=True)
        if best is None or len(ff)<best[0]: best=(len(ff),ns,sorted(av))
        if not ff: break
    return best

if __name__=='__main__':
    a=int(sys.argv[1]); b=int(sys.argv[2])
    r=solve_pair(a,b)
    if r:
        n,ns,av=r
        v=E.forward(ns)
        out=f'full_{a}_{b}_{39033-n}.json'
        json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open(out,'w'))
        json.dump({str(k):str(int(x)) for k,x in ns.items()}, open(out.replace('.json','_seed.json'),'w'))
        print("WROTE",out,"score",39033-n)
