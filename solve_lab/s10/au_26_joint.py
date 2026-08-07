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
w=list(v0); fwd2(w,8)
avw=L.all_atom_values(w)
def nz(v):
    av=L.all_atom_values(v); return [a for a in range(L.NA) if av[a] and a not in SS], set(L.failing_eqs(av)), av
print('x_24548 in atoms:', sorted(L.var_atoms[24548]), flush=True)
for a in sorted(L.var_atoms[24548]): print(f'   a{a}: neq={len(L.atom2eq.get(a,{}))} out={L.atom_out.get(a)} src={L.atom_src[a][:90]}')
print('x_14853 in atoms:', sorted(L.var_atoms[14853]))
for a in sorted(L.var_atoms[14853]): print(f'   a{a}: neq={len(L.atom2eq.get(a,{}))} out={L.atom_out.get(a)} src={L.atom_src[a][:90]}')
print(flush=True)
print('=== congruence 2: x_28730 += d together with x_24548 tracking x_25442 ===', flush=True)
for d in (1, 7, 12345):
    w2=list(w); w2[28730]+=d; fwd2(w2)
    dd = w2[25442]-w[25442]
    w2[24548]+=dd; fwd2(w2)
    nzo,f,av=nz(w2)
    print(f'  d={d}: dx25442={dd}  nonzero-outside-seven={nzo}  score={L.NEQ-len(f)} out12={len(f-E)} {sorted(f-E)[:10]}', flush=True)
    print(f'     A1 mod p moved: {(av[22230]-avw[22230])%P!=0}', flush=True)
print('=== congruence 1: x_7068 += d together with x_14853 tracking x_1308 ===', flush=True)
for d in (1, 7, 12345):
    w2=list(w); w2[7068]+=d; fwd2(w2)
    dd = w2[1308]-w[1308]
    w2[14853]+=dd; fwd2(w2)
    nzo,f,av=nz(w2)
    C0n=(av[22229]+7376877*av[35762])%P
    print(f'  d={d}: dx1308={dd}  nonzero-outside-seven={nzo}  score={L.NEQ-len(f)} out12={len(f-E)} {sorted(f-E)[:10]}', flush=True)
    print(f'     C0 moved: {C0n != (avw[22229]+7376877*avw[35762])%P}', flush=True)
