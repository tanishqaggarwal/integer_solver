import sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast, intsolve
from flint import fmpz_mat
D=pickle.load(open(sys.argv[1],'rb'))
r0=D['r0']; cols=D['cols']; nonlin=D['nonlin']; base=D['base']; rounds=D['rounds']
lo=int(sys.argv[2]); hi=int(sys.argv[3])
out=sys.argv[4] if len(sys.argv)>4 else 'js3_out.json'
for maxr in range(lo,hi+1):
    S=[]
    for r in sorted(rounds):
        if r<=maxr: S+=rounds[r]
    Sset=set(S)
    atoms=set(r0)
    for f in S: atoms|=set(cols[f])
    nl={a for f,a in nonlin if f in Sset}
    use=sorted(a for a in atoms if a not in nl)
    print(f"[round<={maxr}] vars={len(S)} atoms={len(atoms)} linear_rows={len(use)}",flush=True)
    A=[[cols[f].get(a,0) for f in S] for a in use]
    b=[-r0.get(a,0) for a in use]
    t0=time.time()
    M=fmpz_mat(A); Mb=fmpz_mat([A[i]+[b[i]] for i in range(len(use))])
    ra=M.rank(); rb=Mb.rank()
    print(f"   rank(A)={ra} rank([A|b])={rb} ({time.time()-t0:.0f}s) rational-feasible={ra==rb}",flush=True)
    if ra!=rb: continue
    t0=time.time()
    sol,ker=intsolve.solve_int(A,b)
    print(f"   HNF {time.time()-t0:.0f}s integer-feasible={sol is not None} kernel={len(ker)}",flush=True)
    if sol is None: continue
    v0=E.forward(base)
    ns=dict(base)
    nmv=0
    for f,d in zip(S,sol):
        if d: ns[f]=v0[f]+d; nmv+=1
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"   EXACT: moved {nmv} vars -> fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:15]}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open(out,'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open(out.replace('.json','_seed.json'),'w'))
    break
