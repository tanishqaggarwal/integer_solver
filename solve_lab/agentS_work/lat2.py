"""Expose the exact obstruction, widen the knob set, and re-solve at every boolean image point."""
import sys, json, collections, pickle, time, math
sys.path.insert(0,'.')
import common as C, lattice as L
import harness as H, engine as E, fast, sparse, intsolve
P=C.P

def system(seed, extra_knobs=()):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    cone=sorted(set(C.cluster_cone())|set(C.CLUSTERKN)|set(extra_knobs))
    T=L.delta_table(v0,bad0,cone)
    aff={f:d for f,(d,a) in T.items() if a and d and not C.isbool(f)}
    atoms=set(bad0)
    for d in aff.values(): atoms|=set(d)
    cache={}
    hs=L.handles_for(v0,bad0,sorted(atoms),cache)
    for a,(f,s) in hs.items():
        if f not in aff: aff[f]={a:s}
    # widen: free vars in the cones of no-handle atoms
    nh=[a for a in atoms if a not in hs]
    add=set()
    for a in nh:
        try: add|=set(E.cone(a)[1])
        except Exception: pass
    add-=set(aff)|set(cone)
    if add:
        T2=L.delta_table(v0,bad0,sorted(add))
        for f,(d,afq) in T2.items():
            if afq and d and not C.isbool(f): aff[f]=d
        atoms=set(bad0)
        for d in aff.values(): atoms|=set(d)
        hs2=L.handles_for(v0,bad0,sorted(atoms-set(hs)),cache)
        hs.update(hs2)
        for a,(f,s) in hs2.items():
            if f not in aff: aff[f]={a:s}
        atoms=set(bad0)
        for d in aff.values(): atoms|=set(d)
    atoms=sorted(atoms)
    return v0,bad0,aff,atoms,hs

def solve(seed,label,extra=()):
    v0,bad0,aff,atoms,hs=system(seed,extra)
    nh=[a for a in atoms if a not in hs]
    knobs=sorted(aff)
    rows=[{f:aff[f][a] for f in knobs if a in aff[f]} for a in atoms]
    rhs=[-bad0.get(a,0) for a in atoms]
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=atoms,verbose=False,maxcore=600,maxcorebits=400_000)
    print(f"[{label}] bad={sorted(bad0)} knobs={len(aff)} atoms={len(atoms)} no-handle={nh} -> {msg}",flush=True)
    return v0,bad0,aff,atoms,hs,sol,msg,rows,rhs,knobs

if __name__=='__main__':
    seed=dict(C.BASE)
    v0,bad0,aff,atoms,hs,sol,msg,rows,rhs,knobs=solve(seed,'cfg0-wide')
    if sol is None:
        # find the maximal feasible subset of rows, and identify the blocking rows
        keep=[]; block=[]
        for i,a in enumerate(atoms):
            idx=keep+[i]
            s2,_,_=sparse.solve_sparse([rows[j] for j in idx],[rhs[j] for j in idx],verbose=False,maxcore=600,maxcorebits=400_000)
            if s2 is not None: keep=idx
            else: block.append(a)
        print("max feasible rows: %d/%d ; BLOCKING atoms: %s"%(len(keep),len(atoms),block))
        s3,_,_=sparse.solve_sparse([rows[j] for j in keep],[rhs[j] for j in keep],verbose=False,maxcore=600,maxcorebits=400_000)
        if s3 is not None:
            ns=dict(seed)
            for f,dv in s3.items():
                if dv: ns[f]=v0[f]+dv
            v=E.forward(ns); bad=E.badatoms(v); ff=E.eqfails(bad)
            print("partial-applied: bad=%s fails=%d SCORE=%d"%(sorted(bad),len(ff),39033-len(ff)))
            json.dump({str(k):str(int(x)) for k,x in ns.items()},open('lat2_seed.json','w'))
            if len(ff)<28:
                json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('S_lat2_%d.json'%(39033-len(ff)),'w'))
                print("wrote S_lat2_%d.json"%(39033-len(ff)))
