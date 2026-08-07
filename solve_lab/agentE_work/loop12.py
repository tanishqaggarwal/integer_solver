"""Monotone chain repair: never re-touch a variable already set (no reverting)."""
import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
MUXV={22162,30213}
def mux(s,n=4):
    v=E.forward(s)
    for _ in range(n):
        s[22162]=v[13682]; s[30213]=v[18956]-v[32237]; v=E.forward(s)
    return v
def score11(s):
    s=dict(s); v=mux(s); av=E.badatoms(v); f=E.eqfails(av)
    return len(f),av,v,s
def run(abit,bbit,maxsteps=80,log=sys.stdout):
    cur={18956:C,abit:1,bbit:1}
    frozen={18956,abit,bbit}|MUXV
    n,av,v,cur=score11(cur)
    best=(n,dict(cur),dict(av))
    print(f"start fails={n} bad={sorted(av)}",file=log,flush=True)
    for step in range(maxsteps):
        if not av: break
        cands=[]
        for aid in av:
            order,fr,seen=E.cone(aid)
            if len(fr)>500: continue
            for f in fr:
                if f in frozen: continue
                val=E.solve_for(aid,f,cur)
                if val is None or val==cur.get(f,0): continue
                cands.append((aid,f,val))
        bm=None
        for aid,f,val in cands:
            s=dict(cur); s[f]=val
            nn,aav,vv,ss=score11(s)
            kv=(len(aav),nn)
            if bm is None or kv<bm[0]: bm=(kv,ss,aav,nn,aid,f,val)
        if bm is None:
            print(f"step{step}: STUCK bad={sorted(av)}",file=log,flush=True); break
        cur=bm[1]; av=bm[2]; n=bm[3]; frozen.add(bm[5])
        print(f"step{step}: a{bm[4]} x_{bm[5]}<-{str(bm[6])[:20]} bad={len(av)} fails={n} {sorted(av)}",file=log,flush=True)
        if n<best[0]:
            best=(n,dict(cur),dict(av))
            vv=E.forward(cur)
            json.dump({f"x_{i}":vv[i] for i in range(E.NV) if vv[i]!=0}, open(f'chain_{abit}_{bbit}_best.json','w'))
    return best
if __name__=='__main__':
    a=int(sys.argv[1]); b=int(sys.argv[2])
    best=run(a,b,int(sys.argv[3]) if len(sys.argv)>3 else 80)
    print("BEST fails",best[0],"bad",sorted(best[2]))
    v=E.forward(best[1])
    json.dump({f"x_{i}":v[i] for i in range(E.NV) if v[i]!=0}, open(f'chain_{a}_{b}_{best[0]}.json','w'))
