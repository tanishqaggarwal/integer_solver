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
def cone(u):
    seen=set(); st=[u]
    while st:
        t=st.pop()
        if t in seen: continue
        seen.add(t)
        a=L.definer.get(t)
        if a is None: continue
        for x in L.avars[a]:
            if x!=t and x not in seen: st.append(x)
    return seen
c=cone(19964)
print('x_19964 cone:',len(c),sorted(c))
fr=sorted(u for u in c if u not in L.definer)
print('free inputs in cone:',fr)
print('definer of x_19964: a1461', L.atom_src[1461])
# price each free input in the cone (raw, frame2)
for u in fr:
    for d in (1,):
        w2=list(w); w2[u]+=d; fwd2(w2)
        av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
        print(f'  x_{u} (+{d}, val={str(w[u])[:24]}): score={L.NEQ-len(f)} out12={len(f-E)}  d19964={(w2[19964]-w[19964])%P!=0} d19964any={w2[19964]!=w[19964]}')
