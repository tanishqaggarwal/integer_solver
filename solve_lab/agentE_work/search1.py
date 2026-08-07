import sys, json, time, pickle, heapq, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E

def key(seed): return tuple(sorted(seed.items()))

def expand(seed, av, limit=None):
    """Generate candidate seeds: solve each bad atom for each free var in its cone."""
    cands=[]
    for aid in sorted(av, key=lambda a: len(E.cone(a)[1])):
        order,fr,seen=E.cone(aid)
        for f in fr:
            val=E.solve_for(aid,f,seed)
            if val is None: continue
            if val==seed.get(f,0): continue
            s=dict(seed); s[f]=val
            cands.append((aid,f,val,s))
    return cands

def run(start, rounds=8, beam=4, log=sys.stdout):
    n0,av0,v0=E.score(start)
    frontier=[(len(av0),n0,start,av0)]
    seen={key(start)}
    best=(n0,start,av0)
    for rd in range(rounds):
        newf=[]
        for nb,nf,seed,av in frontier:
            t0=time.time()
            cands=expand(seed,av)
            print(f"[rd{rd}] from bad={sorted(av)} fails={nf}: {len(cands)} cands ({time.time()-t0:.0f}s)",file=log,flush=True)
            scored=[]
            for aid,f,val,s in cands:
                k=key(s)
                if k in seen: continue
                nn,aav,vv=E.score(s)
                scored.append((len(aav),nn,s,aav,aid,f,val))
            scored.sort(key=lambda x:(x[0],x[1]))
            for row in scored[:beam*3]:
                nbb,nnf,s,aav,aid,f,val=row
                print(f"   cand a{aid} x_{f}<- {str(val)[:25]} => bad={nbb} fails={nnf} {sorted(aav)}",file=log,flush=True)
                seen.add(key(s)); newf.append((nbb,nnf,s,aav))
                if nnf<best[0]: best=(nnf,s,aav)
        if not newf: break
        newf.sort(key=lambda x:(x[0],x[1]))
        frontier=newf[:beam]
        print(f"[rd{rd}] best so far fails={best[0]} bad={sorted(best[2])}",file=log,flush=True)
        if best[0]==0: break
    return best

if __name__=='__main__':
    C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    start={} if len(sys.argv)<2 else json.load(open(sys.argv[1]))
    start={int(k):int(v) for k,v in start.items()} if start else {}
    best=run(start, rounds=6, beam=3)
    print("BEST fails",best[0],"bad",sorted(best[2]))
    v=E.forward(best[1])
    json.dump({f"x_{i}":v[i] for i in range(E.NV) if v[i]!=0}, open('search1_best.json','w'))
    json.dump({str(k):str(vv) for k,vv in best[1].items()}, open('search1_seed.json','w'))
