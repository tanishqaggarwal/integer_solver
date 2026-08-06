"""S10 step 82: is the residual GLOBALLY pinned, or only locally?

Every rigidity result so far -- the closed 79-column system, the certificates, the
sacrifice analysis -- is a LINEARISATION AT ONE POINT.  Linear algebra at a point
says nothing about far-away points of a nonlinear system.

Test it the blunt way: randomise the non-boolean free inputs (in blocks and
wholesale), forward-evaluate, and see whether the failing CHECK set changes at all.
If the same six checks fail from every starting point, they are structurally
pinned.  If some randomisation changes the set, there is global freedom that every
local method of mine has been blind to.
"""
import os, sys, random, collections, json, time
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

base = L.load(os.path.join(HERE, 'forward_state.json'))
NONBOOL = [u for u in ad.FREE if u not in BOOL]
print(f'non-boolean free inputs: {len(NONBOOL)}')
av0 = L.all_atom_values(base)
F0 = [a for a in range(L.NA) if av0[a] and a not in atom_out]
f0 = L.failing_eqs(av0)
print(f'base: failing checks {F0}, failing equations {len(f0)}')

rng = random.Random(12345)
seen = collections.Counter()
best = (len(f0), None, None)
t0 = time.time()
TRIALS = [('1 random input', 1), ('10 random inputs', 10),
          ('100 random inputs', 100), ('1000 random inputs', 1000),
          ('ALL non-boolean', len(NONBOOL))]
for tag, k in TRIALS:
    for rep in range(4):
        v = list(base)
        picks = rng.sample(NONBOOL, k) if k < len(NONBOOL) else NONBOOL
        for u in picks:
            v[u] = rng.randrange(1, P)
        ad.fwd(v, rounds=3)
        av = L.all_atom_values(v)
        checks = tuple(a for a in range(L.NA) if av[a] and a not in atom_out)
        fail = L.failing_eqs(av)
        seen[checks] += 1
        if len(fail) < best[0]:
            best = (len(fail), tag, checks)
        if rep == 0:
            print(f'  {tag:<20} -> {len(checks)} failing checks, '
                  f'{len(fail)} failing equations   checks={list(checks)[:10]}',
                  flush=True)
print(f'\ndistinct failing-check sets seen: {len(seen)}  ({time.time()-t0:.0f}s)')
for cks, n in seen.most_common(6):
    print(f'  x{n}: {len(cks)} checks {list(cks)[:12]}')
print(f'\nbest failing-equation count over all randomisations: {best[0]} '
      f'(base {len(f0)})')

# Does the six-check core ALWAYS survive?
core = {7930, 29539, 35759, 35760}
always = all(core <= set(c) for c in seen)
print(f'\ndo 7930, 29539, 35759, 35760 fail in EVERY randomisation? {always}')
if not always:
    for c in seen:
        if not core <= set(c):
            print(f'  *** a randomisation avoided part of the core: {list(c)[:14]}')
            break
