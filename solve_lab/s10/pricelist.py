"""S10 step 57: the PRICE LIST, and an equation-scored beam search.

Two mistakes in my earlier greedy searches:
  * they scored by NUMBER OF NONZERO ATOMS, not failing equations;
  * they refused any move that increased the atom count, even when the new atoms
    were 1-equation checks and the closed one cost 15.

Correct objective: minimise failing EQUATIONS, and treat every check as a purchase
with a price = how many equations it lives in.  Budget to beat: 7.
"""
import os, sys, collections, json, time, heapq
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
base = L.load(os.path.join(HERE, 'forward_state.json'))

price = {a: len(L.atom2eq.get(a, {})) for a in range(L.NA) if a not in atom_out}
hist = collections.Counter(price.values())
print('price distribution over all 10,792 checks (equations per check):')
for k in sorted(hist)[:12]:
    print(f'   {k:>3} equations : {hist[k]:>5} checks')
cheap = [a for a, c in price.items() if c <= 2]
print(f'checks costing <= 2 equations: {len(cheap)}')

av0 = L.all_atom_values(base)
FAIL0 = [a for a in range(L.NA) if av0[a] and a not in atom_out]
print(f'\nforward-frame failing checks and their prices:')
for a in FAIL0:
    print(f'   a{a:<7} price={price[a]}')
print(f'   total equation slots {sum(price[a] for a in FAIL0)}, '
      f'actual failing {len(L.failing_eqs(av0))}')


def fwd(v, rounds=3):
    for _ in range(rounds):
        for x in ad.ORDER:
            nv = T.solve_lin(L.definer[x], x, v)
            if nv is not None:
                v[x] = nv
    return v


def score(v):
    av = L.all_atom_values(v)
    return len(L.failing_eqs(av)), av


# ---- beam search, scored by FAILING EQUATIONS -------------------------------
BEAM = int(os.environ.get('BEAM', 12))
DEPTH = int(os.environ.get('DEPTH', 8))
f0, av0 = score(base)
frontier = [(f0, base, ())]
best = (f0, base, ())
seen = set()
t0 = time.time()
for depth in range(DEPTH):
    nxt = []
    for f, v, path in frontier:
        av = L.all_atom_values(v)
        checks = [a for a in range(L.NA) if av[a] and a not in atom_out]
        for a in checks:
            for u in sorted(L.avars[a]):
                if u in L.definer or (a, u) in path:
                    continue
                nv = T.solve_lin(a, u, v)
                if nv is None or nv == v[u]:
                    continue
                w = list(v); w[u] = nv; fwd(w)
                fw, aw = score(w)
                key = tuple(sorted(x for x in range(L.NA) if aw[x]))
                if key in seen:
                    continue
                seen.add(key)
                nxt.append((fw, w, path + ((a, u),)))
                if fw < best[0]:
                    best = (fw, w, path + ((a, u),))
                    print(f'  depth {depth}: NEW BEST failing={fw} '
                          f'score={L.NEQ-fw} via a{a}/x_{u}', flush=True)
    nxt.sort(key=lambda t: t[0])
    frontier = nxt[:BEAM]
    if not frontier:
        break
    print(f'depth {depth}: {len(nxt)} branches, best={frontier[0][0]} '
          f'({time.time()-t0:.0f}s)', flush=True)

print(f'\nBEST failing={best[0]} score={L.NEQ-best[0]}  path={best[2]}')
if best[0] < f0:
    T.save(best[1], os.path.join(HERE, f'price_{L.NEQ-best[0]}.json'))
    print(f'saved price_{L.NEQ-best[0]}.json')
