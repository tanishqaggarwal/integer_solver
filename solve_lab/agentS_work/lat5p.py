"""Full exact integer solve at every BFS image configuration.  Re-measures knobs, handles AND
   targets at each configuration (all three are configuration-dependent)."""
import os
import sys, json, collections, pickle, time
sys.path.insert(0,'.')
import common as C, lat2
import harness as H, engine as E, fast, sparse
P=C.P
seen=pickle.load(open('bfs_image.pkl','rb'))
print("configurations:",len(seen),flush=True)
best=(10**9,None,None)
WK=int(os.environ["WK"]); NW=int(os.environ["NW"])
for i,(k,a) in enumerate(sorted(seen.items(), key=lambda kv:len(kv[1]))):
    if i<18 or i%NW!=WK: continue
    seed=dict(C.BASE); seed.update(a)
    t0=time.time()
    try:
        v0,bad0,aff,atoms,hs=lat2.system(seed)
    except Exception as e:
        print("cfg%d: system ERR %s"%(i,type(e).__name__),flush=True); continue
    nh=[x for x in atoms if x not in hs]
    knobs=sorted(aff)
    rows=[{f:aff[f][x] for f in knobs if x in aff[f]} for x in atoms]
    rhs=[-bad0.get(x,0) for x in atoms]
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=atoms,verbose=False,maxcore=600,maxcorebits=400_000)
    tag='FEASIBLE!!!' if sol is not None else msg
    print("cfg%d (|on|=%d) bad=%s knobs=%d atoms=%d nohandle=%s -> %s (%.0fs)"%(
          i,len(a),sorted(bad0),len(knobs),len(atoms),nh,tag,time.time()-t0),flush=True)
    if sol is not None:
        ns=dict(seed)
        for f,dv in sol.items():
            if dv: ns[f]=v0[f]+dv
        v=E.forward(ns); bad=E.badatoms(v); ff=E.eqfails(bad)
        print("   APPLIED: bad=%s fails=%d SCORE=%d"%(sorted(bad),len(ff),39033-len(ff)),flush=True)
        if len(ff)<best[0]:
            best=(len(ff),dict(ns),v)
            json.dump({"x_%d"%j:int(v[j]) for j in range(E.NV) if v[j]!=0},open('S_lat5_%d.json'%(39033-len(ff)),'w'))
            print("   WROTE S_lat5_%d.json"%(39033-len(ff)),flush=True)
print("BEST",39033-best[0] if best[1] else None,flush=True)
