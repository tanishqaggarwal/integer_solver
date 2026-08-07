"""Joint feasibility for a SUBSET of activated bits (pin system only, CORE selector atoms excluded)."""
import sys, time, pickle, itertools, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B, sparse
C=B.C
def feas_subset(bits, maxr=5, maxv=4000):
    base={18956:C}
    for b in bits: base[b]=1
    frozen={18956}|set(bits)
    v0,bad0,S,cols,nonlin,rounds=B.build_closure(base,frozen,maxr=maxr,maxv=maxv)
    for mr in sorted(rounds):
        Sp=[]
        for r in sorted(rounds):
            if r<=mr: Sp+=rounds[r]
        use,rows,rhs,nd=B.system(bad0,Sp,cols,nonlin)
        sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=False)
        if sol is not None:
            ns=dict(base)
            for f,d in sol.items():
                if d: ns[f]=v0[f]+d
            v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
            return True,msg,len(Sp),len(use),(len(ff),sorted(set(av)-B.CORE)),ns
    return False,msg,len(Sp),len(use),None,None
if __name__=='__main__':
    R=pickle.load(open('scan_B.pkl','rb'))
    fe=sorted(b for b,d in R.items() if d['feas'])
    rnd=random.Random(7)
    print("feasible B bits:",len(fe),flush=True)
    print("--- PAIRS ---",flush=True)
    pairs=[tuple(rnd.sample(fe,2)) for _ in range(12)]
    for p in pairs:
        t0=time.time()
        ok,msg,nv,nr,ex,_=feas_subset(list(p))
        print(f"pair {p}: FEAS={ok} vars={nv} rows={nr} exact={ex} msg={msg[:60]} ({time.time()-t0:.0f}s)",flush=True)
    print("--- TRIPLES ---",flush=True)
    for _ in range(6):
        t=tuple(rnd.sample(fe,3)); t0=time.time()
        ok,msg,nv,nr,ex,_=feas_subset(list(t))
        print(f"triple {t}: FEAS={ok} vars={nv} rows={nr} exact={ex} msg={msg[:60]} ({time.time()-t0:.0f}s)",flush=True)
