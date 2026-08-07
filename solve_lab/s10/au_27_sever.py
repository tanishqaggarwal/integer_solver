import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]; SS=set(SEVEN)
E12=frozenset([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
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

# ---- (A) compensator analysis for an arbitrary atom X: atoms whose eq-set is inside eqs(X) ----
def compensators(X):
    EX=frozenset(L.atom2eq[X])
    return [a for a in range(L.NA) if a!=X and frozenset(L.atom2eq[a])<=EX]
for X in (3576, 21617, 29539, 7930, 37887, 40826, 41512, 37662):
    c=compensators(X)
    print(f'a{X}: neq={len(L.atom2eq[X])} out12={len(set(L.atom2eq[X])-E12)} compensators-inside={c}')
print()

# ---- (B) sever points ----
# which variables does x_7068 change, and which of them feed a29539 ?
def changed_by(u, d=12345):
    w2=list(w); w2[u]+=d; fwd2(w2)
    return [t for t in range(L.NVARS) if w2[t]!=w[t]]
ch7068=set(changed_by(7068)); ch28730=set(changed_by(28730))
print('vars changed by x_7068:',len(ch7068))
print('vars changed by x_28730:',len(ch28730), sorted(ch28730))

def path_to(target_vars, changed):
    """variables in `changed` that are ancestors of any var in target_vars"""
    anc=set(); st=list(target_vars)
    while st:
        t=st.pop()
        if t in anc: continue
        anc.add(t)
        a=L.definer.get(t)
        if a is None: continue
        for x in L.avars[a]:
            if x not in anc: st.append(x)
    return sorted(anc & changed)
p1 = path_to(L.avars[29539], ch7068)
p2 = path_to(L.avars[7930], ch28730)
print('\nsever candidates on x_7068 -> a29539 :', p1)
for t in p1:
    a=L.definer.get(t)
    if a is None: print(f'   x_{t}: FREE'); continue
    print(f'   x_{t}: definer a{a} neq={len(L.atom2eq.get(a,{}))} out12={len(set(L.atom2eq.get(a,{}))-E12)} src={L.atom_src[a][:70]}')
print('\nsever candidates on x_28730 -> a7930 :', p2)
for t in p2:
    a=L.definer.get(t)
    if a is None: print(f'   x_{t}: FREE'); continue
    print(f'   x_{t}: definer a{a} neq={len(L.atom2eq.get(a,{}))} out12={len(set(L.atom2eq.get(a,{}))-E12)} src={L.atom_src[a][:70]}')
