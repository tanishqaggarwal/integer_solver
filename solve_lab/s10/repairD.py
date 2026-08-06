"""S10 step 17: shift x_7068 by k*p and explicitly close the two atoms the ripple
leaves broken (29539, 40826) using their handles."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
base = L.load(BEST)

for a in (29539, 40826):
    print(f'\na{a}: {L.atom_src[a][:400]}')
    print('   vars:', sorted(L.avars[a]))
    for u in sorted(L.avars[a]):
        free = u not in L.definer
        print(f'     x_{u:<7} val={str(base[u])[:26]:<28} free={free} '
              f'natoms={len(L.var_atoms[u])} neqs={len(L.var_eqs[u])}')

print('\n=== shift test with explicit handle repair ===')
for k in (1, 2, 7, -3228258):
    v = list(base)
    L.ripple(v, {7068: base[7068] + k * P}, block=BLOCK)
    av = L.all_atom_values(v)
    broken = [a for a in (29539, 40826) if av[a] != 0]
    fixed = []
    for a in broken:
        opts = []
        for u in sorted(L.avars[a]):
            nv = T.solve_lin(a, u, v)
            if nv is None or nv == v[u]:
                continue
            opts.append((len(L.var_eqs[u]), u, nv))
        opts.sort()
        if opts:
            _, u, nv = opts[0]
            L.ripple(v, {u: nv}, block=BLOCK)
            fixed.append((a, u))
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    D = v[7068] - v[2099]
    print(f'k={k:<10} repaired={fixed} nz={nz} failing={len(fail)} '
          f'score={L.NEQ-len(fail)} D%7376877={D%7376877}')
