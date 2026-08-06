"""S10 step 60: the TRUE cost of sacrificing a hub check.

Checks 2423, 31670 and 19297 each hit FOUR of the six inconsistency certificates,
so a hitting set costs  cost(X) + 1 + 1  (certs 0 and 5 have cost-1 members).
Budget 7  =>  we win if some hub X costs <= 4 after compensation.

"Price" (number of equations the atom lives in) is only an UPPER bound: the
39,026 witness saves 5 of its 12 equations by cancellation.  Compute, per hub,
how many of its equations can be saved by atoms in them that carry an adjustable
handle and whose own footprint stays inside the sacrificed region.
"""
import os, sys, collections, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
v = L.load(os.path.join(HERE, 'forward_state.json'))
av = L.all_atom_values(v)
FAIL = set(a for a in range(L.NA) if av[a] and a not in atom_out)
FAILQ = set(L.failing_eqs(av))

# adjustable atoms: those carrying a solo free handle, with its granularity
SOLO = collections.defaultdict(list)
for u in range(L.NVARS):
    if u not in L.definer and len(L.var_atoms[u]) == 1:
        SOLO[L.var_atoms[u][0]].append(u)
GRAN = {}
for a, us in SOLO.items():
    g = 0
    for u in us:
        r = T.lin_parts(a, u, v)
        if r:
            g = math.gcd(g, abs(r[0]))
    if g:
        GRAN[a] = g
print(f'atoms with an adjustable handle: {len(GRAN)}')

HUBS = [2423, 31670, 19297, 21617, 31672, 33929, 34397, 40065]
print(f'{"hub":>8} {"price":>6} {"eqs already failing":>20} {"compensable":>12} '
      f'{"TRUE COST":>10}')
for X in HUBS:
    E = sorted(L.atom2eq.get(X, ()))
    if not E:
        continue
    already = [e for e in E if e in FAILQ]
    saveable = 0
    detail = []
    for e in E:
        if e in FAILQ:
            continue
        m, sq, co = L.eq_atoms[e]
        # atoms in e (other than X) that are adjustable AND whose own equations
        # lie inside E (so using them costs nothing new)
        helpers = [a for a in co
                   if a != X and a in GRAN
                   and set(L.atom2eq.get(a, ())) <= set(E) | FAILQ]
        loose = [a for a in co if a != X and a in GRAN]
        if helpers:
            saveable += 1
            detail.append((e, helpers[:3]))
    print(f'{X:>8} {len(E):>6} {len(already):>20} {saveable:>12} '
          f'{len(E)-len(already)-saveable:>10}')
    for e, h in detail[:4]:
        print(f'          eq {e}: helpers {h}')

print('\n=== how many adjustable atoms sit in each hub equation at all? ===')
for X in (2423, 31670, 19297):
    E = sorted(L.atom2eq.get(X, ()))
    print(f'\nhub a{X} ({len(E)} equations):')
    for e in E:
        m, sq, co = L.eq_atoms[e]
        adj = [a for a in co if a != X and a in GRAN]
        nz = [a for a in co if av[a]]
        print(f'   eq {e:<6} atoms={len(co):<3} adjustable={len(adj):<3} '
              f'nonzero_now={nz} failing={e in FAILQ}')
