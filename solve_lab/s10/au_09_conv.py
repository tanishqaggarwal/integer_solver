import os, sys, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8)
t0=time.time(); av=L.all_atom_values(w); print('all_atom_values %.2fs'%(time.time()-t0))
t0=time.time(); f=L.failing_eqs(av); print('failing_eqs %.2fs'%(time.time()-t0))
# convergence test: perturb a var, run with increasing rounds
import random
random.seed(1)
FREE=[u for u in range(L.NVARS) if u not in definer]
print('n free params in frame2:', len(FREE))
for u in [9118, 8731, 6418, 7068, 28730, 642, 29854, 31864]:
    prev=None
    line=[]
    for r in (1,2,3,4,6,8,12,20):
        w2=list(w); w2[u]+=1; fwd2(w2,r)
        av2=L.all_atom_values(w2); s=L.NEQ-len(L.failing_eqs(av2))
        line.append((r,s))
    print(f'  x_{u}+1 score by rounds: {line}')
# also: does fwd2 on the base state with many rounds stay put?
w3=list(v0); fwd2(w3,40); print('base after 40 rounds identical:', w3==v0)
# ORDER: how many vars are in the cyclic (non-topo) part?
print('len(ad.ORDER)',len(ad.ORDER),'topo',len(L.topo),'cyc',len(ad.ORDER)-len(L.topo))
