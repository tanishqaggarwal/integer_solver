import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
print('a29090 (definer of x_2099):', L.atom_src[29090])
print('  eqs of a29090:', len(L.atom2eq.get(29090,{})))
print('  atom_out:', L.atom_out.get(29090))
print()
# ancestor cone of x_2099 restricted to definers
def cone(u, maxn=200000):
    seen=set(); st=[u]
    while st:
        t=st.pop()
        if t in seen: continue
        seen.add(t)
        a=L.definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w!=t and w not in seen: st.append(w)
    return seen
c2099 = cone(2099)
print('x_2099 ancestor cone size:', len(c2099))
free2099 = sorted(u for u in c2099 if u not in L.definer)
print('free inputs in it:', len(free2099), free2099[:40])
for u in free2099[:40]:
    print(f'   x_{u}: val={str(v[u])[:40]} atoms={sorted(L.var_atoms[u])[:8]} n_atoms={len(L.var_atoms[u])}')
