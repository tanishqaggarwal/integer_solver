"""S10 step 85: verify the message sweep's best state by actual construction.

The sweep says the ENTIRE message space (31,008 reachable states) zeroes at most
2 of the 6 failing checks -- 35759 and 35760 -- via a single bit.  The sweep uses
the linear bit model; confirm it by building the state and forward-evaluating.
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
spec = json.load(open(os.path.join(HERE, 'msgsweep.json')))
groups = spec['groups']; cur = spec['cur']; best = spec['best']
CHECKS = spec['checks']
print(f'groups (size, currently set): {[(len(us), c) for (s, us), c in zip(groups, cur)]}')
print(f'best group counts {best} vs current {cur}')

v = L.load(os.path.join(HERE, 'forward_state.json'))
av0 = L.all_atom_values(v)
f0 = L.failing_eqs(av0)
print(f'base: failing {len(f0)}, score {L.NEQ-len(f0)}')

# realise the target counts: turn bits on/off within each group
w = list(v)
changed = []
for gi, ((s, us), c0, c1) in enumerate(zip(groups, cur, best)):
    on = [u for u in us if v[u] == 1]
    off = [u for u in us if v[u] == 0]
    need = c1 - c0
    if need > 0:
        for u in off[:need]:
            w[u] = 1; changed.append((u, 1))
    elif need < 0:
        for u in on[:-need]:
            w[u] = 0; changed.append((u, 0))
print(f'bits changed: {changed}')
ad.fwd(w, rounds=3)
av = L.all_atom_values(w)
nz = [a for a in range(L.NA) if av[a]]
nzc = [a for a in nz if a not in atom_out]
fail = L.failing_eqs(av)
print(f'\nafter: nonzero atoms {len(nz)} (checks {len(nzc)}) failing {len(fail)} '
      f'score {L.NEQ-len(fail)}')
print(f'  failing checks now: {nzc[:20]}')
for c in CHECKS:
    print(f'    a{c}: {"ZEROED" if av[c] == 0 else "still nonzero"}')
if len(fail) < len(f0):
    T.save(w, os.path.join(HERE, f'msg_{L.NEQ-len(fail)}.json'))
    print(f'  saved msg_{L.NEQ-len(fail)}.json')

# and the honest question: does the linear bit model match reality?
print('\nlinear model predicted 35759 and 35760 zeroed; '
      f'reality: 35759 {"zero" if av[35759] == 0 else "NONZERO"}, '
      f'35760 {"zero" if av[35760] == 0 else "NONZERO"}')
