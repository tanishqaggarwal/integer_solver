import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B, fast, sparse
C=B.C
bit=int(sys.argv[1]); MC=int(sys.argv[2]) if len(sys.argv)>2 else 1200
s={18956:C, bit:1}
v=E.forward(s)
for _ in range(6): s[14853]=v[13682]; v=E.forward(s)
v0=E.forward(s); bad0=E.badatoms(v0)
print('bad',sorted(bad0),flush=True)
S=[];cols={};nonlin=set();processed=set();rounds={}
pending=set(bad0)
for rnd in range(6):
    newS=set()
    for a in pending: newS|=set(E.cone(a)[1])
    newS-=set(S)|{18956,bit}
    newS=sorted(newS)
    if not newS: break
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
    processed|=pending; pending=(aff|set(bad0))-processed
    if len(S)>4000: break
for mr in sorted(rounds):
    Sp=[]
    for r in sorted(rounds):
        if r<=mr: Sp+=rounds[r]
    Sset=set(Sp); atoms=set(bad0)
    for f in Sp: atoms|=set(cols[f])
    nl={a for f,a in nonlin if f in Sset}
    use=sorted(a for a in atoms if a not in nl)
    rowmap={a:{} for a in use}
    for f in Sp:
        for a,c in cols[f].items():
            if a in rowmap: rowmap[a][f]=c
    rows=[rowmap[a] for a in use]; rhs=[-bad0.get(a,0) for a in use]
    t0=time.time()
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=use,verbose=True,maxcore=MC)
    print(f'  r<={mr}: vars={len(Sp)} rows={len(use)} -> {msg[:120]} ({time.time()-t0:.0f}s)',flush=True)
    if sol is None: continue
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    vv=E.forward(ns); av=E.badatoms(vv); ff=E.eqfails(av)
    print(f'  EXACT fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:12]}',flush=True)
    json.dump({f"x_{i}":int(vv[i]) for i in range(E.NV) if vv[i]!=0}, open(f'big01_{bit}_{39033-len(ff)}.json','w'))
    if not ff: break
