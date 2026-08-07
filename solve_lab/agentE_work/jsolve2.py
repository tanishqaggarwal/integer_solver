import sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast, intsolve
D=pickle.load(open(sys.argv[1],'rb'))
r0=D['r0']; cols=D['cols']; nonlin=D['nonlin']; base=D['base']; rounds=D['rounds']
maxr=int(sys.argv[2]) if len(sys.argv)>2 else 1
S=[]
for r in sorted(rounds):
    if r<=maxr: S+=rounds[r]
atoms=set(r0)
for f in S: atoms|=set(cols[f])
atoms=sorted(atoms)
nl={a for f,a in nonlin if f in set(S)}
print("vars",len(S),"atoms",len(atoms),"nonlinear atoms",len(nl))
use=[a for a in atoms if a not in nl]
print("using",len(use),"linear rows; dropping",len(atoms)-len(use))
A=[[cols[f].get(a,0) for f in S] for a in use]
b=[-r0.get(a,0) for a in use]
t0=time.time()
sol,ker=intsolve.solve_int(A,b)
print("HNF solve %.1fs"%(time.time()-t0),"feasible:",sol is not None,"kernel",len(ker))
if sol is None:
    keep=[]
    for i in range(len(use)):
        s2,_=intsolve.solve_int([A[j] for j in keep+[i]],[b[j] for j in keep+[i]])
        if s2 is not None: keep.append(i)
    drop=[use[i] for i in range(len(use)) if i not in keep]
    print("greedy max solvable rows",len(keep),"/",len(use)," dropped:",drop[:20])
    sol,ker=intsolve.solve_int([A[j] for j in keep],[b[j] for j in keep])
if sol is not None:
    v0=E.forward(base)
    ch={f:v0[f]+d for f,d in zip(S,sol) if d}
    print("moving",len(ch),"vars")
    ns=dict(base)
    for f,val in ch.items(): ns[f]=val
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print("EXACT after move: fails",len(ff),"score",39033-len(ff),"bad atoms",sorted(av)[:20])
    json.dump({f"x_{i}":v[i] for i in range(E.NV) if v[i]!=0}, open(sys.argv[3] if len(sys.argv)>3 else 'js2_out.json','w'))
    json.dump({str(k):str(x) for k,x in ns.items()}, open('js2_seed.json','w'))
