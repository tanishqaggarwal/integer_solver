"""Minimum-equation-cost sacrifice set: solve the exact integer atom system with a small set of
   atoms DROPPED (allowed to carry residual), preferring atoms with tiny equation footprint."""
import sys, json, collections, pickle, time, itertools
sys.path.insert(0,'.')
import common as C, lat2
import harness as H, engine as E, fast, sparse
P=C.P
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}

def try_cfg(seed, tag, maxdrop=3, verbose=True):
    v0,bad0,aff,atoms,hs=lat2.system(seed)
    knobs=sorted(aff)
    rowd={a:{f:aff[f][a] for f in knobs if a in aff[f]} for a in atoms}
    rhs={a:-bad0.get(a,0) for a in atoms}
    def solve(keep):
        return sparse.solve_sparse([rowd[a] for a in keep],[rhs[a] for a in keep],
                                   verbose=False,maxcore=600,maxcorebits=400_000)[0]
    if solve(atoms) is not None:
        print(f"[{tag}] FULL SYSTEM FEASIBLE",flush=True); return ('FEASIBLE',atoms,None)
    # candidate drop atoms: order by equation footprint ascending
    cand=sorted(atoms, key=lambda a:(NF.get(a,99),a))
    best=None
    for k in range(1,maxdrop+1):
        for S in itertools.combinations(cand[:14],k):
            keep=[a for a in atoms if a not in S]
            sol=solve(keep)
            if sol is None: continue
            ns=dict(seed)
            for f,dv in sol.items():
                if dv: ns[f]=v0[f]+dv
            v=E.forward(ns); bad=E.badatoms(v); ff=E.eqfails(bad)
            cost=len(ff)
            if best is None or cost<best[0]:
                best=(cost,S,dict(ns),v,sorted(bad))
                print(f"[{tag}] drop={S} nf={[NF.get(a) for a in S]} -> bad={sorted(bad)} fails={cost} SCORE={39033-cost}",flush=True)
                if cost<7:
                    json.dump({"x_%d"%j:int(v[j]) for j in range(E.NV) if v[j]!=0},open('S_sac_%d.json'%(39033-cost),'w'))
                    print(f"[{tag}] *** WROTE S_sac_{39033-cost}.json ***",flush=True)
        if best and best[0]<=7: break
    return ('DROP',best,None)

if __name__=='__main__':
    seed=dict(C.BASE)
    print("candidate drop atoms at cfg0 by footprint:")
    v0,bad0,aff,atoms,hs=lat2.system(seed)
    print(sorted([(NF.get(a,99),a) for a in atoms])[:16])
    try_cfg(seed,'cfg0',maxdrop=2)
