"""S10 step 1: exact equation-space picture of the residual at the 39,026 witness.

Reuses the s9 model (validated: 0 mismatches over all 39,033 equations).
"""
import os, sys, json, collections

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L   # chdir's into s9 and loads all caches

BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')

v = L.load(BEST)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a] != 0]
fail = L.failing_eqs(av)
print(f'nonzero atoms : {len(nz)} -> {nz}')
for a in nz:
    print(f'   atom {a}: value={av[a]}')
    print(f'      src = {L.atom_src[a][:200]}')
print(f'failing eqs   : {len(fail)} -> {fail}')
print()

# Which equations does each nonzero atom appear in?
for a in nz:
    eqs = sorted(L.atom2eq.get(a, {}))
    print(f'atom {a} appears in {len(eqs)} equations: {eqs}')
print()

print('=== per failing equation: full atom decomposition ===')
for i in fail:
    m, sq, co = L.eq_atoms[i]
    tot = sum(c * av[a] for a, c in co.items())
    print(f'\neq {i}:  mult={m}  square={sq}  n_atoms={len(co)}  combo={tot}')
    for a, c in sorted(co.items()):
        mark = '  <== NONZERO' if av[a] else ''
        print(f'    coeff {c:>14}  atom {a:>6}  val={av[a]}{mark}')
        if not av[a]:
            print(f'         src={L.atom_src[a][:150]}')
