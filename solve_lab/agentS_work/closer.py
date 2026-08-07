"""Per-state closer: forward -> exact affine lattice solve on the max-feasible row set -> eq score."""
import sys, json, collections, pickle, time
sys.path.insert(0,'.')
import common as C, lat2, lattice as L
import harness as H, engine as E, fast, sparse
P=C.P
def close(seed, drop_try=True, verbose=False):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    if not bad0: return 0,seed,{},v0
    try: v0,bad0,aff,atoms,hs=lat2.system(seed)
    except Exception as e: return None
    knobs=sorted(aff)
    rows=[{f:aff[f][a] for f in knobs if a in aff[f]} for a in atoms]
    rhs=[-bad0.get(a,0) for a in atoms]
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=atoms,verbose=False,maxcore=600,maxcorebits=400_000)
    if sol is None and drop_try:
        # greedily keep the largest feasible prefix set
        keep=[]
        for i in range(len(atoms)):
            idx=keep+[i]
            s2,_,_=sparse.solve_sparse([rows[j] for j in idx],[rhs[j] for j in idx],verbose=False,maxcore=600,maxcorebits=400_000)
            if s2 is not None: keep=idx
        if keep:
            sol,_,_=sparse.solve_sparse([rows[j] for j in keep],[rhs[j] for j in keep],verbose=False,maxcore=600,maxcorebits=400_000)
    ns=dict(seed)
    if sol:
        for f,dv in sol.items():
            if dv: ns[f]=v0[f]+dv
    v=E.forward(ns); bad=E.badatoms(v); ff=E.eqfails(bad)
    return len(ff),ns,bad,v
if __name__=='__main__':
    import glob
    seeds=[('cfg0',{})]
    for fn in sorted(glob.glob('bfs_hit_*.json')):
        seeds.append((fn,{int(k):int(v) for k,v in json.load(open(fn)).items()}))
    best=(10**9,None)
    for tag,ex in seeds:
        s=dict(C.BASE); s.update(ex)
        t0=time.time()
        r=close(s)
        if r is None: print(f"{tag}: ERR"); continue
        n,ns,bad,v=r
        print(f"{tag}: fails={n} SCORE={39033-n} bad={sorted(bad)} ({time.time()-t0:.0f}s)",flush=True)
        if n<best[0]:
            best=(n,ns)
            if n<7:
                json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0},open('S_best_%d.json'%(39033-n),'w'))
                print("  WROTE S_best_%d.json"%(39033-n))
    print("BEST",39033-best[0])
