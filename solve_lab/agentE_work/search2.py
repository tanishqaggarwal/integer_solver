"""Per-bit activation cost: set one OR-tree bit to 1, greedily repair its collateral atoms."""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E

C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
CORE={20212,20215,24403,747}

def cone_of_var(u0):
    seen=set(); st=[u0]; fr=set()
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        dv=E.definer[u]
        if dv is None: fr.add(u); continue
        for w in E.avars[dv[0]]:
            if w!=u: st.append(w)
    return seen,fr

def greedy(seed, frozen, maxsteps=5, log=None):
    n,av,v=E.score(seed)
    cur=dict(seed); curav=dict(av); curn=n
    best=(len(set(curav)-CORE), curn, dict(cur), dict(curav))
    visited={tuple(sorted(cur.items()))}
    for step in range(maxsteps):
        tgts=[a for a in curav if a not in CORE]
        if not tgts: break
        cands=[]
        for aid in tgts:
            order,fr,seen=E.cone(aid)
            if len(fr)>400: continue
            for f in fr:
                if f in frozen: continue
                val=E.solve_for(aid,f,cur)
                if val is None or val==cur.get(f,0): continue
                cands.append((f,val))
        bm=None
        for f,val in cands:
            s=dict(cur); s[f]=val
            k=tuple(sorted(s.items()))
            if k in visited: continue
            nn,aav,vv=E.score(s)
            coll=len(set(aav)-CORE)
            keyv=(coll,nn)
            if bm is None or keyv<bm[0]: bm=(keyv,s,aav,nn)
        if bm is None: break
        cur=bm[1]; curav=bm[2]; curn=bm[3]
        visited.add(tuple(sorted(cur.items())))
        if (len(set(curav)-CORE),curn) < (best[0],best[1]):
            best=(len(set(curav)-CORE),curn,dict(cur),dict(curav))
        if log: print("    step%d coll=%d fails=%d bad=%s"%(step,len(set(curav)-CORE),curn,sorted(curav)),file=log,flush=True)
        if not set(curav)-CORE: break
    return best

if __name__=='__main__':
    roots=[int(x) for x in sys.argv[1:]] or [7715,34554]
    out={}
    log=sys.stdout
    for root in roots:
        seen,fr=cone_of_var(root)
        bits=sorted(fr)
        print("root %d: %d free candidates"%(root,len(bits)),file=log,flush=True)
        for f in bits:
            s={18956:C, f:1}
            n0,av0,_=E.score(s)
            t0=time.time()
            b=greedy(s,{f,18956},maxsteps=5)
            out[(root,f)]=(b[0],b[1],dict(b[2]),sorted(b[3]))
            print("  bit x_%d: start fails=%d bad=%s -> coll=%d fails=%d bad=%s (%.0fs)"%(f,n0,sorted(av0),b[0],b[1],sorted(b[3]),time.time()-t0),file=log,flush=True)
            pickle.dump(out,open('search2.pkl','wb'))
