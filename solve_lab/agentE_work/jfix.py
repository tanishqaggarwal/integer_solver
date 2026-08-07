"""Load closure, solve r0 + J d = 0 over Z (linear atoms only), apply, verify."""
import sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, jsolve as J, intsolve
D=pickle.load(open(sys.argv[1] if len(sys.argv)>1 else 'jacC.pkl','rb'))
r0=D['r0']; cols=D['cols']; nonlin=D['nonlin']; S=D['S']; base=D['base']
badpairs={}
for f,a in nonlin: badpairs.setdefault(a,set()).add(f)
allat=set(r0)
for f in cols: allat|=set(cols[f])
lin_at=[a for a in sorted(allat) if a not in badpairs]
print("vars",len(S),"atoms",len(allat),"linear atoms",len(lin_at))
A=[[cols[f].get(a,0) for f in S] for a in lin_at]
b=[-r0.get(a,0) for a in lin_at]
t0=time.time()
sol,ker=intsolve.solve_int(A,b)
print("solve time %.1fs"%(time.time()-t0),"solution?",sol is not None,"kernel dim",len(ker))
if sol is None:
    # find a maximal solvable subset greedily
    print("infeasible on full linear set; probing subsets")
    keep=[]
    for i,a in enumerate(lin_at):
        cand=keep+[i]
        s2,_=intsolve.solve_int([A[j] for j in cand],[b[j] for j in cand])
        if s2 is not None: keep=cand
    print("max greedy solvable rows:",len(keep),"of",len(lin_at))
    print("dropped atoms:",[lin_at[i] for i in range(len(lin_at)) if i not in keep])
    sol,ker=intsolve.solve_int([A[j] for j in keep],[b[j] for j in keep])
if sol is not None:
    ns=dict(base)
    for f,d in zip(S,sol):
        if d: ns[f]=ns.get(f,0)+d
    v=J.mux(ns); av=E.badatoms(v); f=E.eqfails(av)
    print("APPLIED: fails",len(f),"bad atoms",sorted(av))
    json.dump({f"x_{i}":v[i] for i in range(E.NV) if v[i]!=0}, open('jfix_out.json','w'))
    json.dump({str(k):str(x) for k,x in ns.items()}, open('jfix_seed.json','w'))
    print("wrote jfix_out.json")
