"""bl_eqs: the failing equations -- can another atom CANCEL a residual atom?"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, FORBID
P = 2**256-2**32-977

for name, F, path in (('F2/w39026', F2, os.path.join(LAB,'best','new_instance_partial_39026.json')),
                      ('CANON/mod9118_0', CANON, os.path.join(HERE,'mod9118_0.json'))):
    v = L.load(path); F.fwd(v, rounds=8)
    av = L.all_atom_values(v)
    fail = L.failing_eqs(av)
    nz = [a for a in range(L.NA) if av[a]]
    print(f'\n===== {name}: {len(fail)} failing equations, nonzero atoms {nz} =====')
    for i in fail:
        m, sq, co = L.eq_atoms[i]
        others = sorted(a for a in co if a not in nz)
        print(f'  eq{i}: mult={m} sq={sq} atoms={ {a: co[a] for a in sorted(co)} }')
        print(f'        value={str(L.eq_value(i, av))[:50]}...  other atoms in eq: {len(others)} {others[:14]}')
        # can any other atom in this equation be made nonzero freely?
        for a in others[:14]:
            free_in = sorted(set(L.avars[a]) & F.FREE)
            print(f'          a{a} (coef {co[a]}) out={L.atom_out.get(a)} freevars={free_in[:8]} src={L.atom_src[a][:70]}')
