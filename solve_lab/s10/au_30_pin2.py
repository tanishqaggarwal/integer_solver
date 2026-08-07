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
# congruence-2 pin: 9367949*(x_24548 - x_25442) = x_7927 = x_15616*x_11052 = p*x_11052
# x_25442 = x_22342 + x_10861 ; x_22342 = x_26732 + x_1627 ; x_26732 = x_11425*x_4432
for u in (1627, 10861, 11425, 15616, 11052, 7927, 24548, 8557, 1237):
    d=L.definer.get(u)
    print(f'x_{u}: free={d is None} definer=a{d} val={str(w[u])[:40]} ==p:{w[u]==P} natoms={len(L.var_atoms[u])}'
          + (f'  defsrc={L.atom_src[d][:70]}' if d is not None else ''))
print()
# ancestor cones -> free inputs, and their raw price
def cone(u):
    seen=set(); st=[u]
    while st:
        t=st.pop()
        if t in seen: continue
        seen.add(t)
        a=definer.get(t)
        if a is None: continue
        for x in L.avars[a]:
            if x!=t and x not in seen: st.append(x)
    return seen
targets={'x_1627':1627,'x_10861':10861,'x_11425':11425,'x_24548':24548}
cand=set()
for k,t in targets.items():
    c=cone(t); fr=sorted(u for u in c if u not in definer)
    print(f'{k}: cone {len(c)}, free inputs {fr}')
    cand |= set(fr)
print('\nraw price of each candidate (frame 2):')
avb=L.all_atom_values(w)
for u in sorted(cand):
    w2=list(w); w2[u]+=1000003; fwd2(w2)
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    nzo=[a for a in range(L.NA) if av[a] and a not in SS]
    print(f'  x_{u}: score={L.NEQ-len(f)} out12={len(f-E)} nonzero-outside-seven={nzo[:8]}', flush=True)
