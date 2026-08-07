"""Is the 'cost 0' claim for x9118/x8731 real?  Move each candidate knob in the WITNESS
frame (frame2: x7068,x28730,x29854,x31864,x642 detached), forward-solve, and count the
equations that break outside the twelve."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
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
print('witness in its own frame after fwd: score %d (failing %d)'%(L.NEQ-len(fe0),len(fe0)))
E=set(e for a in env.SEVEN for e in L.atom2eq[a])
SEVEN=env.SEVEN
CANDS=[9118,8731,642,29854,31864,1329,10903,9413,17325,2099,7068,28730,4432,19964,
       950,6947,9629,15120,23754,33168,35619,1613,1844,21574,29305,2892,6090,28355]
print('%-8s %-8s %-8s %s'%('var','cost','score','delta on the seven'))
for u in CANDS:
    v=list(w); v[u]=v[u]+1; fwd(v)
    av=L.all_atom_values(v); fe=set(L.failing_eqs(av))
    cost=len(fe-E)
    d=[av[a]-av0[a] for a in SEVEN]
    ds=', '.join('a%d%+d'%(SEVEN[i],x) if abs(x)<10**9 else 'a%d~%dd'%(SEVEN[i],len(str(abs(x)))) for i,x in enumerate(d) if x)
    print('x%-7d %-8d %-8d %s'%(u,cost,L.NEQ-len(fe),ds))
