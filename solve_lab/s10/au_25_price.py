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
FREESET=set(u for u in range(L.NVARS) if u not in definer)
def fwd2(v,rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8); avb=L.all_atom_values(w)
CHECKS=[a for a in range(L.NA) if a not in {definer[t] for t in definer}]
def nzchecks(v):
    av=L.all_atom_values(v)
    return [a for a in range(L.NA) if av[a] and a not in SS]
print('a7930 vars:', sorted(L.avars[7930]), ' x_24548 free:',24548 in FREESET, ' x_25442 def a', L.definer.get(25442))
print('a29539 vars:', sorted(L.avars[29539]),' x_14853 free:',14853 in FREESET,' x_1308 def a', L.definer.get(1308))
print('x_25442 src:', L.atom_src[L.definer[25442]])
print('x_1308 src:', L.atom_src[L.definer[1308]])
print()
for u,d in [(14853,1),(24548,1),(14853,3),(24548,5)]:
    w2=list(w); w2[u]+=d; fwd2(w2)
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    print(f'x_{u}+{d}: score={L.NEQ-len(f)} out12={len(f-E)} nonzero-outside-seven={nzchecks(w2)}')
    print(f'   out eqs: {sorted(f-E)}')
