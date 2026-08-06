"""S10 step 6: trace the ripple from x_7068 exactly -- which vars move, by how much,
and which atoms end up broken.  Then test the compensated move (move the mirror too)."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
CLEAN = os.path.join(HERE, 'clean_state.json')
v0 = L.load(CLEAN)
DELTA = v0[2099] - v0[7068]        # what we need x_7068 to move by
print('required delta =', DELTA)
print('delta mod p    =', DELTA % P)

def run(seeds, block=(), tag=''):
    v = L.load(CLEAN)
    changed, _ = L.ripple(v, dict(seeds), block=block)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    return v, changed, nz, fail

v, changed, nz, fail = run({7068: v0[7068] + DELTA})
print(f'\nripple from x_7068: {len(changed)} vars changed, {len(nz)} nz atoms, '
      f'{len(fail)} failing')
print('\nvariables that moved by exactly DELTA (1:1 carriers):')
one2one = [u for u in changed if v[u] - v0[u] == DELTA]
print('  ', sorted(one2one))
print('variables that moved by -DELTA:')
print('  ', sorted([u for u in changed if v[u] - v0[u] == -DELTA]))
print('variables that moved by other amounts:')
other = {u: v[u]-v0[u] for u in changed if abs(v[u]-v0[u]) != abs(DELTA)}
for u, d in list(other.items())[:25]:
    r = d / DELTA if DELTA else 0
    print(f'   x_{u:<7} d={str(d)[:34]:<36} d/DELTA={r if abs(r)<1e12 else "big"}')
print(f'   ... {len(other)} total')

print('\nbroken atoms and the handle available in each:')
for a in nz:
    print(f'  a{a}: {L.atom_src[a][:150]}')
    for m, c in sorted(L.polys[a].items()):
        if len(m) == 2:
            vals = [v[w] for w in m]
            tag = ''
            if P in vals: tag = ' <== has p factor (handle: quotient is free)'
            if vals == [0, 0]: tag = ' <== DORMANT'
            print(f'      monomial x_{m[0]}*x_{m[1]}  vals=({str(vals[0])[:18]},{str(vals[1])[:18]}){tag}')
