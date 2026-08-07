"""Greedy repair with the (1,1)-branch MUX fixpoint maintained."""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E

C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
MUXV={22162,30213}

def mux(s, n=4):
    v=E.forward(s)
    for _ in range(n):
        s[22162]=v[13682]; s[30213]=v[18956]-v[32237]
        v=E.forward(s)
    return v

def score11(s):
    s=dict(s); v=mux(s)
    av=E.badatoms(v); f=E.eqfails(av)
    return len(f), av, v, s

def greedy(seed, frozen, maxsteps=25, log=sys.stdout, tag=''):
    n,av,v,cur=score11(seed)
    best=(len(av), n, dict(cur), dict(av))
    visited={tuple(sorted(cur.items()))}
    print(f"{tag} start fails={n} bad={sorted(av)}",file=log,flush=True)
    for step in range(maxsteps):
        cands=[]
        for aid in av:
            order,fr,seen=E.cone(aid)
            if len(fr)>500: continue
            for f in fr:
                if f in frozen or f in MUXV: continue
                val=E.solve_for(aid,f,cur)
                if val is None or val==cur.get(f,0): continue
                cands.append((aid,f,val))
        bm=None
        for aid,f,val in cands:
            s=dict(cur); s[f]=val
            k=tuple(sorted(s.items()))
            if k in visited: continue
            nn,aav,vv,ss=score11(s)
            kv=(len(aav),nn)
            if bm is None or kv<bm[0]: bm=(kv,ss,aav,nn,aid,f,val)
        if bm is None:
            print(f"{tag} step{step}: no move",file=log,flush=True); break
        cur=bm[1]; av=bm[2]; n=bm[3]
        visited.add(tuple(sorted(cur.items())))
        print(f"{tag} step{step}: a{bm[4]} x_{bm[5]}<-{str(bm[6])[:22]} => bad={len(av)} fails={n} {sorted(av)}",file=log,flush=True)
        if (len(av),n)<(best[0],best[1]): best=(len(av),n,dict(cur),dict(av))
        if not av: break
    return best

if __name__=='__main__':
    abit=int(sys.argv[1]); bbit=int(sys.argv[2])
    seed={18956:C, abit:1, bbit:1}
    b=greedy(seed, {abit,bbit,18956}, maxsteps=int(sys.argv[3]) if len(sys.argv)>3 else 25, tag=f'[{abit},{bbit}]')
    print("BEST bad=",b[0],"fails=",b[1],"badset",sorted(b[3]))
    v=E.forward(b[2])
    out=f'loop11_{abit}_{bbit}_{b[1]}.json'
    json.dump({f"x_{i}":v[i] for i in range(E.NV) if v[i]!=0}, open(out,'w'))
    json.dump({str(k):str(x) for k,x in b[2].items()}, open(f'seed_{abit}_{bbit}.json','w'))
    print("wrote",out)
