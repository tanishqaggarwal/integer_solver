import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
for u in (28599, 22665, 28961, 17499, 7075):
    a=L.definer.get(u)
    print(f'x_{u}: definer a{a}  src={L.atom_src[a] if a is not None else None}')
    print(f'   n_eqs of a{a} = {len(L.atom2eq.get(a,{}))}, atom_out={L.atom_out.get(a)}')
    # who else uses x_u
    print(f'   x_{u} appears in atoms: {sorted(L.var_atoms[u])[:20]} (n={len(L.var_atoms[u])})')
    print()
# how many variables equal exactly p?
eqp=[u for u in range(L.NVARS) if v[u]==P]
print('variables equal to p exactly:', len(eqp))
print('  sample', eqp[:30])
