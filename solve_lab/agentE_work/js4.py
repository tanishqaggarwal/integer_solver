import sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, sparse
D=pickle.load(open(sys.argv[1],'rb'))
r0=D['r0']; cols=D['cols']; nonlin=D['nonlin']; base=D['base']; rounds=D['rounds']
lo=int(sys.argv[2]); hi=int(sys.argv[3]); out=sys.argv[4] if len(sys.argv)>4 else 'js4.json'
for maxr in range(lo,hi+1):
    S=[]
    for r in sorted(rounds):
        if r<=maxr: S+=rounds[r]
    Sset=set(S); atoms=set(r0)
    for f in S: atoms|=set(cols[f])
    nl={a for f,a in nonlin if f in Sset}
    use=sorted(a for a in atoms if a not in nl)
    print(f"[<= {maxr}] vars={len(S)} atoms={len(atoms)} linrows={len(use)}",flush=True)
    rowmap={a:{} for a in use}
    for f in S:
        for a,c in cols[f].items():
            if a in rowmap: rowmap[a][f]=c
    rows=[rowmap[a] for a in use]; rhs=[-r0.get(a,0) for a in use]
    t0=time.time()
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=use)
    print(f"  -> {msg} ({time.time()-t0:.1f}s)",flush=True)
    if sol is None: continue
    v0=E.forward(base); ns=dict(base); nmv=0
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d; nmv+=1
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"  EXACT: moved {nmv} -> fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:15]}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open(out,'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open(out.replace('.json','_seed.json'),'w'))
    break
