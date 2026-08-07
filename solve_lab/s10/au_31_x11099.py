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
KEY=[1627,10861,11425,25442,22342,26732,24548,19964,2099,7068,28730,4432,14853,1308,6418,12553]
for u,d in [(11099,1),(11099,1000003)]:
    w2=list(w); w2[u]+=d; fwd2(w2)
    ch=[t for t in range(L.NVARS) if w2[t]!=w[t]]
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    print(f'x_{u} += {d}: vars changed = {len(ch)} {ch[:20]}  score={L.NEQ-len(f)}')
    print('   key vars moved mod p:', {k:( (w2[k]-w[k])%P!=0 ) for k in KEY if w2[k]!=w[k]})
print()
print('x_11099 appears in atoms:')
for a in sorted(L.var_atoms[11099]):
    print(f'   a{a}: out={L.atom_out.get(a)} neq={len(L.atom2eq.get(a,{}))} src={L.atom_src[a][:110]}')
print('x_34660:', w[34660], 'definer', L.definer.get(34660))
print('a10937 src:', L.atom_src[10937])
print('x_37413:', w[37413], 'definer', L.definer.get(37413))
