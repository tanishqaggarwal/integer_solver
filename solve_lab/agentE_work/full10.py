"""Joint (1,0) solve: activate ONE a-tree bit, close the MUX with x_12186 / x_16742,
   then solve the whole remaining atom system as one exact linear Diophantine system."""
import sys, time, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, full11 as F
C=F.C
def mux10(s,n=6):
    v=E.forward(s)
    for _ in range(n):
        s[12186]=v[13682]; s[16742]=v[18956]-v[32237]; v=E.forward(s)
    return v,s
def solve_bit(abit, verbose=True, maxr=6, iters=4):
    s={18956:C, abit:1}
    v0,s=mux10(s)
    best=None
    for it in range(iters):
        bad0=E.badatoms(v0)
        if verbose: print(f"[{abit}] iter{it} bad={sorted(bad0)} fails={len(E.eqfails(bad0))}",flush=True)
        if not bad0:
            return 0,s,[]
        S,cols,nonlin,rounds=F.closure(v0,bad0,{18956,abit},maxr=maxr,verbose=verbose)
        got=None
        for mr in sorted(rounds):
            Sp=[]
            for r in sorted(rounds):
                if r<=mr: Sp+=rounds[r]
            Sset=set(Sp); atoms=set(bad0)
            for f in Sp: atoms|=set(cols[f])
            nl={a for f,a in nonlin if f in Sset}
            use=sorted(a for a in atoms if a not in nl)
            rowmap={a:{} for a in use}
            for f in Sp:
                for a,c in cols[f].items():
                    if a in rowmap: rowmap[a][f]=c
            rows=[rowmap[a] for a in use]; rhs=[-bad0.get(a,0) for a in use]
            import sparse
            sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=False)
            if verbose: print(f"   r<={mr}: vars={len(Sp)} rows={len(use)} nl={len(nl)} -> {msg[:70]}",flush=True)
            if sol is None: continue
            ns=dict(s)
            for f,d in sol.items():
                if d: ns[f]=v0[f]+d
            v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
            if verbose: print(f"   EXACT fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:12]}",flush=True)
            if got is None or len(ff)<got[0]: got=(len(ff),ns,sorted(av),v)
            if not ff: break
        if got is None:
            return best
        if best is None or got[0]<best[0]: best=(got[0],dict(got[1]),got[2])
        if got[0]==0: return best
        s=dict(got[1]); v0=got[3]
    return best
if __name__=='__main__':
    a=int(sys.argv[1])
    r=solve_bit(a)
    if r:
        n,ns,av=r
        v=E.forward(ns)
        out=f'f10_{a}_{39033-n}.json'
        json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open(out,'w'))
        json.dump({str(k):str(int(x)) for k,x in ns.items()}, open(out.replace('.json','_seed.json'),'w'))
        print("WROTE",out,"score",39033-n,"bad",av)
