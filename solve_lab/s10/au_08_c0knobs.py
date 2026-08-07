import os, sys, collections, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = set([2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])
DETACH = {7068:22229, 28730:22230, 29854:35758, 31864:35761, 642:35762}
definer = {t:a for t,a in L.definer.items() if t not in DETACH}
ORDER = [t for t in ad.ORDER if t not in DETACH]
def fwd2(v, rounds=8):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None: v[u]=nv
    return v
v0 = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w = list(v0); fwd2(w)
av0 = L.all_atom_values(w)
t0=time.time(); fwd2(list(w)); print('one fwd2 (8 rounds) takes %.1fs'%(time.time()-t0))
base2099 = w[2099]%P
base7068 = w[7068]%P
base28730 = w[28730]%P

for u,d in [(2081,1),(4287,1),(6418,1),(31861,1),(9118,1),(2081,7),(4287,3),(6418,999),(31861,5)]:
    w2=list(w); w2[u]+=d; fwd2(w2)
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    out = sorted(f-E)
    print(f'x_{u} += {d}: score={L.NEQ-len(f)}  failing_outside_12={len(out)} {out[:12]}'
          f'  dx2099modp={(w2[2099]-w[2099])%P!=0}  dx7068modp={(w2[7068]-w[7068])%P!=0}'
          f'  dx28730modp={(w2[28730]-w[28730])%P!=0}')
