import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

def rep(v, tag):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    f = L.failing_eqs(av)
    print(f'{tag}: score {L.NEQ-len(f)}  nonzero {len(nz)} {nz[:12]}')
    return av, nz, f

t0=time.time()
w26 = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
rep(w26, 'witness39026')
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
rep(v0, 'mod9118_0')
print(f'load+eval {time.time()-t0:.1f}s')

# timing of the cheap pre-filter
t=time.time(); v=list(v0); ad.fwd(v, rounds=6); print(f'ad.fwd rounds=6: {time.time()-t:.2f}s')
t=time.time(); av=L.all_atom_values(v); print(f'all_atom_values: {time.time()-t:.2f}s')
t=time.time(); f=L.failing_eqs(av); print(f'failing_eqs: {time.time()-t:.2f}s')
t=time.time(); v=list(v0); ad.fwd(v, rounds=3); print(f'ad.fwd rounds=3: {time.time()-t:.2f}s')

# what does fwd do to the 39026 witness in the canonical frame?
w=list(w26); ad.fwd(w, rounds=6); rep(w,'witness39026 after canonical fwd')

# frame-2 forward
DETACH = {7068:22229, 28730:22230, 29854:35758, 31864:35761, 642:35762}
def2 = {t_:a for t_,a in definer.items() if t_ not in DETACH}
ORDER2 = [t_ for t_ in ad.ORDER if t_ not in DETACH]
def fwd2(v, rounds=6):
    for _ in range(rounds):
        for u in ORDER2:
            nv = T.solve_lin(def2[u], u, v)
            if nv is not None: v[u]=nv
    return v
w=list(w26); fwd2(w); rep(w,'witness39026 after FRAME2 fwd')
t=time.time(); w=list(w26); fwd2(w, rounds=6); print(f'fwd2 rounds=6: {time.time()-t:.2f}s')
