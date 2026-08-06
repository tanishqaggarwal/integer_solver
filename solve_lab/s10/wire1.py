"""S10 step 42: drive the whole p-wire to a different value.

The 220-variable wire is 219 copies of one root, x_26064, held by a single BARE
pin  a37694 = x_26064 - p  which appears in only 12 equations.  Every handle in
the instance enters as  wire * handle, so the wire is the sole reason all 1,249
handles are p-quantised.  Set the wire to w != p (blocking a37694) and the
handles become granularity w -- for w = 1, completely unquantised, and every
pin  bit*(x_B - HUGE) = s*x_C  becomes absorbable.

Cost: a37694 itself, in 12 equations.  Question: what is left after that?
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
ROOT, ROOTATOM = 26064, 37694
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIRE = [u for u in range(L.NVARS) if base[u] == P]


def fwd_block(v, block, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            a = L.definer[u]
            if a in block:
                continue
            nv = T.solve_lin(a, u, v)
            if nv is not None:
                v[u] = nv
    return v


print(f'wire size {len(WIRE)}; root x_{ROOT} pinned by a{ROOTATOM} '
      f'(in {len(L.atom2eq.get(ROOTATOM, {}))} equations)')

for W in (1, 0, -1, 2):
    v = list(base)
    v[ROOT] = W
    for u in WIRE:
        v[u] = W
    fwd_block(v, {ROOTATOM})
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    nzc = [a for a in nz if a not in atom_out]
    nzg = [a for a in nz if a in atom_out]
    fail = L.failing_eqs(av)
    stillwire = sum(1 for u in WIRE if v[u] == W)
    print(f'\nwire = {W}: wire members holding = {stillwire}/{len(WIRE)}')
    print(f'   nonzero atoms {len(nz)} (checks {len(nzc)}, gates {len(nzg)}) '
          f'failing={len(fail)} score={L.NEQ-len(fail)}')
    print(f'   checks: {nzc[:30]}')
    if nzg:
        print(f'   broken gates: {nzg[:15]}')
    if W == 1:
        T.save(v, os.path.join(HERE, 'wire1_state.json'))

# how p-quantised are the handles once the wire is 1?
print('\n=== handle granularity with wire = 1 ===')
v = L.load(os.path.join(HERE, 'wire1_state.json'))
solo = collections.defaultdict(list)
for u in range(L.NVARS):
    if u in L.definer:
        continue
    ats = L.var_atoms[u]
    if len(ats) == 1:
        solo[ats[0]].append(u)
import math
kinds = collections.Counter()
for a, us in solo.items():
    g = 0
    for u in us:
        r = T.lin_parts(a, u, v)
        if r is not None:
            g = math.gcd(g, abs(r[0]))
    if g == 0:
        kinds['rigid'] += 1
    elif g == 1:
        kinds['FREE (granularity 1)'] += 1
    elif g % P == 0:
        kinds['p-quantised'] += 1
    else:
        kinds[f'other'] += 1
for k, n in kinds.most_common():
    print(f'   {k:<28} {n}')
