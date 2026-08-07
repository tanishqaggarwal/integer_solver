"""Exact integer linear system over the affine cone knobs (handles included).
   Rows = atoms in play, cols = knobs, solve D n = -R over Z.  Any solution is then
   VERIFIED by full re-propagation (affinity/superposition is not assumed, it is checked)."""
import sys, json, collections, pickle, time, math
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast, sparse
P=C.P

def delta_table(v0, bad0, knobs):
    T={}
    for f in knobs:
        o=v0[f]
        b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2})
        keys=set(b1)|set(bad0)
        d1={a:b1.get(a,0)-bad0.get(a,0) for a in keys}; d1={a:x for a,x in d1.items() if x}
        aff=all(b2.get(a,0)-bad0.get(a,0)==2*d1.get(a,0) for a in keys)
        T[f]=(d1,aff)
    return T

def handles_for(v0,bad0,atoms,cache):
    """pure single-atom affine knobs for the given atoms (searched in each atom's cone)."""
    out={}
    for a in atoms:
        if a in cache: 
            if cache[a]: out[a]=cache[a]
            continue
        try: order,fr,seen=E.cone(a)
        except Exception: cache[a]=None; continue
        found=None
        for f in fr:
            o=v0[f]
            try:
                b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2})
            except Exception: continue
            keys=set(b1)|set(bad0)
            d1={x:b1.get(x,0)-bad0.get(x,0) for x in keys}; d1={x:y for x,y in d1.items() if y}
            if len(d1)==1 and a in d1 and b2.get(a,0)-bad0.get(a,0)==2*d1[a]:
                if found is None or abs(d1[a])<abs(found[1]): found=(f,d1[a])
        cache[a]=found
        if found: out[a]=found
    return out

def build_and_solve(seed, label, extra_knobs=(), verbose=True):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    cone=sorted(set(C.cluster_cone())|set(C.CLUSTERKN)|set(extra_knobs))
    T=delta_table(v0,bad0,cone)
    aff={f:d for f,(d,a) in T.items() if a and d and not C.isbool(f)}
    atoms=set(bad0)
    for d in aff.values(): atoms|=set(d)
    cache={}
    hs=handles_for(v0,bad0,sorted(atoms),cache)
    # add handles as knobs
    for a,(f,s) in hs.items():
        if f not in aff: aff[f]={a:s}
    atoms=set(bad0)
    for d in aff.values(): atoms|=set(d)
    atoms=sorted(atoms)
    nohandle=[a for a in atoms if a not in hs]
    if verbose:
        print(f"[{label}] bad={sorted(bad0)}  affine knobs={len(aff)}  atoms in play={len(atoms)}  no-handle atoms={nohandle}",flush=True)
    knobs=sorted(aff)
    rows=[{f:aff[f].get(a,0) for f in knobs if a in aff[f]} for a in atoms]
    rhs=[-bad0.get(a,0) for a in atoms]
    t0=time.time()
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=atoms,verbose=verbose,maxcore=600,maxcorebits=200_000)
    print(f"[{label}] solve_sparse -> {msg} ({time.time()-t0:.0f}s)",flush=True)
    return v0,bad0,aff,atoms,hs,sol,msg

if __name__=='__main__':
    seed=dict(C.BASE)
    v0,bad0,aff,atoms,hs,sol,msg=build_and_solve(seed,'cfg0')
    if sol:
        ns=dict(seed)
        for f,dv in sol.items():
            if dv: ns[f]=v0[f]+dv
        v=E.forward(ns); bad=E.badatoms(v); ff=E.eqfails(bad)
        print("APPLIED (full re-propagation): bad=%s fails=%d SCORE=%d"%(sorted(bad),len(ff),39033-len(ff)))
        if len(ff)<7:
            json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('S_lat_%d.json'%(39033-len(ff)),'w'))
            print("WROTE S_lat_%d.json"%(39033-len(ff)))
