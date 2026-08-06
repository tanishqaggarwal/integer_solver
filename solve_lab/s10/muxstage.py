"""S10 step 9: staged construction on the MUX branch, measuring after every stage."""
import os, sys, json, re
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977

def measure(v, tag, show=True):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'[{tag}] nz_atoms={len(nz)} failing={len(fail)} score={L.NEQ-len(fail)}')
    if show:
        for a in nz:
            print(f'    a{a:<6} neq={len(L.atom2eq.get(a,{})):<3} {L.atom_src[a][:100]}')
    return av, nz, fail

v = L.load(os.path.join(HERE, 'mux_on.json'))
measure(v, 'mux_on', show=False)

# ---- stage 1: satisfy the two load pins that x_4287 lit, on FREE inputs -------
print('\n--- stage 1: pin x_31861 and x_14865 ---')
for a, u in ((3568, 31861), (3570, 14865)):
    print(f'  a{a} full: {L.atom_src[a]}')
    nv = L.solve_for(a, u, v)
    print(f'    solve x_{u} -> {str(nv)[:60]}  (free={u not in L.definer})')
    if nv is not None:
        ch, _ = L.ripple(v, {u: nv})
        print(f'    rippled {len(ch)} vars')
measure(v, 'stage1', show=True)
json.dump({f'x_{i}': v[i] for i in range(L.NVARS) if v[i] != 0},
          open(os.path.join(HERE, 'mux_s1.json'), 'w'))
