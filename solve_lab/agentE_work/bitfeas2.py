"""Per-bit pin-system integer feasibility, with round-by-round growth and binding-row report."""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast, sparse
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
CORE={20212,20215,24403,747}
MAXR=int(__import__('os').environ.get('MAXR','5'))
MAXV=int(__import__('os').environ.get('MAXV','3000'))

def cone_free(a):
    return set(E.cone(a)[1])

def build_closure(base, frozen, maxr=MAXR, maxv=MAXV, verbose=False):
    v0=E.forward(base); bad0=E.badatoms(v0)
    S=[]; cols={}; nonlin=set(); processed=set(); rounds={}
    pending=set(bad0)-CORE
    for rnd in range(maxr+1):
        newS=set()
        for a in pending: newS|=cone_free(a)
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
        processed|=pending
        pending=(aff|set(bad0))-processed-CORE
        if verbose: print(f"    rnd{rnd}: +{len(newS)} vars ({len(S)} tot) {time.time()-t0:.0f}s",flush=True)
        if len(S)>maxv: break
    return v0,bad0,S,cols,nonlin,rounds

def system(bad0,S,cols,nonlin):
    Sset=set(S); atoms=set(bad0)-CORE
    for f in S: atoms|=set(cols[f])
    atoms-=CORE
    nl={a for f,a in nonlin if f in Sset}
    use=sorted(a for a in atoms if a not in nl)
    rowmap={a:{} for a in use}
    for f in S:
        for a,c in cols[f].items():
            if a in rowmap: rowmap[a][f]=c
    return use,[rowmap[a] for a in use],[-bad0.get(a,0) for a in use],len(nl)

def analyse(base, frozen, verbose=False):
    v0,bad0,S,cols,nonlin,rounds=build_closure(base,frozen,verbose=verbose)
    out=[]
    for maxr in sorted(rounds):
        Sp=[]
        for r in sorted(rounds):
            if r<=maxr: Sp+=rounds[r]
        use,rows,rhs,ndrop=system(bad0,Sp,cols,nonlin)
        sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=False)
        out.append((maxr,len(Sp),len(use),ndrop,msg,sol))
        if sol is not None: break
    return bad0,out,S,cols,nonlin,rounds,v0

if __name__=='__main__':
    which=sys.argv[1]
    bits=[int(x) for x in sys.argv[2:]]
    res={}
    for b in bits:
        t0=time.time()
        base={18956:C, b:1}
        bad0,out,S,cols,nonlin,rounds,v0=analyse(base,{18956,b})
        maxr,nv,nr,nd,msg,sol=out[-1]
        feas = sol is not None
        nz=sum(1 for x in (sol or {}).values() if x)
        print(f"bit x_{b}: pins={sorted(set(bad0)-CORE)} rounds<= {maxr} vars={nv} rows={nr} dropped_nonlin={nd} "
              f"FEASIBLE={feas} moved={nz} msg={msg[:90]} ({time.time()-t0:.0f}s)",flush=True)
        res[b]=(feas,msg,sorted(set(bad0)-CORE),nv,nr,{k:int(x) for k,x in (sol or {}).items() if x})
        pickle.dump(res,open(f'bitfeas_{which}.pkl','wb'))
