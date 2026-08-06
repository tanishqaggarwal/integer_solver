import os, sys, collections, json, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]

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
f0 = L.failing_eqs(av0)
print('frame2 fwd of witness: score', L.NEQ-len(f0), 'failing', f0, 'identical to v0:', w==v0)

def report(tag, w2):
    av = L.all_atom_values(w2)
    f = L.failing_eqs(av)
    nz = [a for a in range(L.NA) if av[a] and a not in set(SEVEN)]
    print(f'[{tag}] score={L.NEQ-len(f)} failing={sorted(f)}')
    print(f'      nonzero atoms outside the SEVEN: {nz[:20]}')
    A=[av[a] for a in SEVEN]
    C0=(A[0]+7376877*A[6])%P
    print(f'      C0 mod p = {C0}')
    print(f'      A1 mod p = {A[1]%P}')
    print(f'      x_2099 mod p = {w2[2099]%P}   x_7068 mod p = {w2[7068]%P}')
    return av,f

av_b,_ = report('base', w)

for d in (1, 2, 12345):
    w2 = list(w); w2[9118] += d; fwd2(w2)
    print(f'\n--- x_9118 += {d} ---')
    report(f'x9118+{d}', w2)
    print('   delta x_2099 =', (w2[2099]-w[2099]))
    print('   delta x_2099 mod p =', (w2[2099]-w[2099])%P)
