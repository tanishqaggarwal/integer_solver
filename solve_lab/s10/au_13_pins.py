import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]; SS=set(SEVEN)
E=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8); avb=L.all_atom_values(w)
for u,d,tag in [(7068,1,'+1'),(7068,P,'+p'),(7068,3*P,'+3p'),
                (28730,1,'+1'),(28730,P,'+p'),(28730,7*P,'+7p'),
                (17325,1,'+1'),(1329,1,'+1'),(10903,1,'+1'),(9413,1,'+1')]:
    w2=list(w); w2[u]+=d; fwd2(w2)
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    nzout=[a for a in range(L.NA) if av[a] and a not in SS]
    ch=[t for t in range(L.NVARS) if w2[t]!=w[t]]
    print(f'x_{u}{tag}: score={L.NEQ-len(f)} out12={len(f-E)} {sorted(f-E)[:10]}')
    print(f'      nonzero atoms outside seven: {nzout}  vars changed: {len(ch)}')
