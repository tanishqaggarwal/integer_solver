"""S10 step 84: EXHAUSTIVE sweep of the entire message space.

bitgroups.py: only 128 boolean inputs move any failing check, and they carry just
5 distinct signature vectors, with multiplicities 75, 50, 1, 1, 1.  Within a group
the bits are interchangeable, so only the COUNT matters:

    reachable message states = 76 * 51 * 2 * 2 * 2 = 31,008

not 2^256.  That is the whole "256-bit codeword" the earlier sessions treated as a
combinatorial wall, and it is enumerable in a second.  Sweep it and report, for
every state, how many of the failing checks are zeroed.
"""
import os, sys, collections, itertools, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
BOOL = set()
for _a, _poly in enumerate(L.polys):
    _ks = list(_poly.items())
    if len(_ks) == 2:
        _sq = [m for m, c in _ks if len(m) == 2 and m[0] == m[1]]
        _li = [m for m, c in _ks if len(m) == 1]
        if _sq and _li and _sq[0][0] == _li[0][0]:
            BOOL.add(_li[0][0])

v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
BFREE = set(u for u in ad.FREE if u in BOOL)
CHECKS = [a for a in range(L.NA) if av[a] and a not in atom_out]
print(f'failing checks: {CHECKS}')

mats = {}
for c in CHECKS:
    g = ad.grad(c, vm)
    mats[c] = {u: d % P for u, d in g.items() if u in BFREE and d % P}
allbits = sorted(set().union(*[set(m) for m in mats.values()]))
sig = collections.defaultdict(list)
for u in allbits:
    sig[tuple(mats[c].get(u, 0) for c in CHECKS)].append(u)
groups = sorted(sig.items(), key=lambda kv: -len(kv[1]))
print(f'\n{len(allbits)} bits, {len(groups)} signature groups: '
      f'{[len(g[1]) for g in groups]}')
for s, us in groups:
    on = sum(1 for u in us if v[u] == 1)
    print(f'  group of {len(us):>3} (currently {on} set): '
          f'{[str(x)[:14] for x in s]}')

# residual of each check now, and how each group's count shifts it
res = [(-av[c]) % P for c in CHECKS]      # amount each check must move by
cur = [sum(1 for u in us if v[u] == 1) for s, us in groups]
print(f'\ncurrent group counts: {cur}')
print(f'required shifts: {[str(r)[:22] for r in res]}')

ranges = [range(0, len(us) + 1) for s, us in groups]
total = 1
for r in ranges: total *= len(r)
print(f'\nsweeping {total} message states ...')
best = None
hist = collections.Counter()
for combo in itertools.product(*ranges):
    zeroed = 0
    for ci, c in enumerate(CHECKS):
        d = 0
        for gi, (s, us) in enumerate(groups):
            d += (combo[gi] - cur[gi]) * s[ci]
        if (av[c] + d) % P == 0:
            zeroed += 1
    hist[zeroed] += 1
    if best is None or zeroed > best[0]:
        best = (zeroed, combo)
print(f'histogram of checks zeroed: {dict(sorted(hist.items()))}')
print(f'\nBEST: {best[0]} of {len(CHECKS)} checks zeroed at group counts {best[1]}')
print(f'  (current state zeroes {sum(1 for c in CHECKS if av[c] % P == 0)})')
if best[0] > 0:
    zc = []
    for ci, c in enumerate(CHECKS):
        d = 0
        for gi, (s, us) in enumerate(groups):
            d += (best[1][gi] - cur[gi]) * s[ci]
        if (av[c] + d) % P == 0:
            zc.append(c)
    print(f'  checks zeroed: {zc}')
    json.dump({'groups': [[list(map(str, s)), us] for s, us in groups],
               'cur': cur, 'best': list(best[1]), 'zeroed': best[0],
               'checks': CHECKS},
              open(os.path.join(HERE, 'msgsweep.json'), 'w'))
    print('  saved msgsweep.json')
