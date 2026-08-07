"""Frame-2 generator scan over the ancestor cone of x2099 (and x28730's drivers):
which variables move alpha0 = a22229 or alpha1 = a22230 mod p, and at what cost?"""
import sys, json, collections, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s10')
import tools as T, ad
P=env.P
DETACH={7068:22229, 28730:22230, 29854:35758, 31864:35761, 642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd(v,rounds=6):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
base=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
w=list(base); fwd(w)
av0=L.all_atom_values(w); fe0=set(L.failing_eqs(av0))
E=set(e for a in env.SEVEN for e in L.atom2eq[a])
SEVEN=env.SEVEN
print('base %d'%(L.NEQ-len(fe0)),flush=True)
def cone(t):
    seen=set(); st=[t]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        a=definer.get(u)
        if a is None: continue
        for x in L.avars[a]:
            if x!=u: st.append(x)
    return seen
C=cone(2099)|cone(28730)|{37158,25297}
C=sorted(C)
print('cone of x2099 u x28730: %d variables'%len(C),flush=True)
res=[]
t0=time.time()
for k,u in enumerate(C):
    v=list(w); v[u]=v[u]+1; fwd(v)
    av=L.all_atom_values(v); fe=set(L.failing_eqs(av))
    cost=len(fe-E)
    d=[av[a]-av0[a] for a in SEVEN]
    if any(d) or cost==0:
        res.append((cost,u,[str(x) for x in d]))
        ds=', '.join('a%d%+d'%(SEVEN[i],x) if abs(x)<10**9 else 'a%d~%dd'%(SEVEN[i],len(str(abs(x)))) for i,x in enumerate(d) if x)
        if any(d): print('  x%-7d cost=%-4d score=%-6d %s'%(u,cost,L.NEQ-len(fe),ds),flush=True)
    if k%50==0: print('   ...%d/%d [%.0fs]'%(k,len(C),time.time()-t0),flush=True)
json.dump(res,open('/home/user/integer_solver/solve_lab/agentA_work/scan2099.json','w'))
print('done',flush=True)
