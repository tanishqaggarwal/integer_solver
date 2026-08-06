import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
for a in (3576, 3578, 3580, 3582, 40005, 1462, 1466, 1463):
    print(f'a{a}: out={L.atom_out.get(a)} neq={len(L.atom2eq.get(a,{}))}')
    print('  SRC:', L.atom_src[a])
    print('  polys:', dict(list(L.polys[a].items())[:8]))
    print()
for u in (26777, 13458, 19247, 7133, 25607, 2081, 31861, 14865, 32010):
    d=L.definer.get(u)
    print(f'x_{u}: free={d is None} definer=a{d} val={str(v[u])[:44]} n_atoms_using={len(L.var_atoms[u])} atoms={sorted(L.var_atoms[u])[:12]}')
