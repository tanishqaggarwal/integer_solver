"""S10 step 76: exhaustive budget search -- ANY sacrifice set costing <= 6 ?

The give-up cost is 7, so any set of checks whose combined equation count is <= 6
and whose removal makes the closed GF(p) system consistent BEATS the deliverable.
There are ~20 price-1 checks in the closed system, so subsets up to size 6 are
directly enumerable: C(20,6) = 38,760 consistency tests at ~5 ms each.

This settles the last open door with an exhaustive answer rather than a search.
"""
import os, sys, collections, itertools, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad
from newton import BOOL

P = ad.P
atom_out = L.atom_out
v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
FAIL = set(a for a in range(L.NA) if av[a] and a not in atom_out)

gcache = {}
def gradc(c):
    if c not in gcache:
        gcache[c] = {u: d for u, d in ad.grad(c, vm).items() if u not in BOOL}
    return gcache[c]

reach_cache = {}
def reach(u):
    if u in reach_cache: return reach_cache[u]
    seen, frontier, checks = {u}, [u], set()
    while frontier:
        nxt = []
        for w in frontier:
            for a in L.var_atoms[w]:
                if ad.dpart(a, w, vm) == 0: continue
                if a not in atom_out:
                    checks.add(a); continue
                t = atom_out[a][1]
                if t == w or t in seen or ad.dpart(a, t, vm) == 0: continue
                seen.add(t); nxt.append(t)
        frontier = nxt
    reach_cache[u] = checks
    return checks

cols = set()
for c in FAIL: cols |= set(gradc(c))
for _ in range(10):
    rs = set()
    for u in cols: rs |= reach(u)
    new = set(cols)
    for c in rs: new |= set(gradc(c))
    if len(new) == len(cols): break
    cols = new
rs = set()
for u in cols: rs |= reach(u)
COLS = sorted(cols)
ROWS = []
for c in sorted(rs):
    g = {u: d for u, d in gradc(c).items() if u in cols}
    if not g and (av[c] % P) == 0: continue
    ROWS.append(c)
n = len(COLS); cidx = {u: i for i, u in enumerate(COLS)}
FULL = []
for c in ROWS:
    r = [0] * (n + 1)
    for u, d in gradc(c).items():
        if u in cidx: r[cidx[u]] = d % P
    r[n] = (-av[c]) % P
    FULL.append(r)
price = {c: len(L.atom2eq.get(c, {})) for c in ROWS}
eqs = {c: set(L.atom2eq.get(c, ())) for c in ROWS}
print(f'closed system {len(ROWS)} rows x {n} cols')
byprice = collections.Counter(price.values())
print('price histogram:', dict(sorted(byprice.items())[:8]))


def consistent(dropidx):
    M = [FULL[i][:] for i in range(len(FULL)) if i not in dropidx]
    m = len(M); rr = 0
    for col in range(n):
        k = next((i for i in range(rr, m) if M[i][col]), None)
        if k is None: continue
        M[rr], M[k] = M[k], M[rr]
        inv = pow(M[rr][col], -1, P)
        M[rr] = [x * inv % P for x in M[rr]]
        for i in range(m):
            if i != rr and M[i][col]:
                f = M[i][col]; Mr = M[rr]
                M[i] = [(M[i][j] - f * Mr[j]) % P for j in range(n + 1)]
        rr += 1
        if rr == m: break
    return not any(M[i][n] for i in range(rr, m))


cheap = [i for i in range(len(ROWS)) if price[ROWS[i]] <= 3]
cheap.sort(key=lambda i: price[ROWS[i]])
print(f'rows with price <= 3: {len(cheap)}  '
      f'-> {[(ROWS[i], price[ROWS[i]]) for i in cheap[:24]]}')

t0 = time.time()
found = None
tests = 0
for size in range(1, 7):
    for S in itertools.combinations(cheap, size):
        cost = len(set().union(*[eqs[ROWS[i]] for i in S]))
        if cost > 6:
            continue
        tests += 1
        if consistent(set(S)):
            found = (cost, [ROWS[i] for i in S])
            break
    print(f'  size {size}: {tests} sets tested within budget, '
          f'{"FOUND" if found else "none"}  ({time.time()-t0:.0f}s)', flush=True)
    if found:
        break

if found:
    cost, atoms = found
    print(f'\n*** SACRIFICE SET WITHIN BUDGET: {atoms} at cost {cost} '
          f'=> potential score {L.NEQ - cost}')
    json.dump({'atoms': atoms, 'cost': cost},
              open(os.path.join(HERE, 'budget6.json'), 'w'))
else:
    print(f'\nEXHAUSTIVE: no sacrifice set of price<=3 checks costs <= 6 equations '
          f'and restores consistency ({tests} sets tested).')
    # what IS the cheapest consistent set among these?
    best = None
    for size in range(1, 5):
        for S in itertools.combinations(cheap, size):
            cost = len(set().union(*[eqs[ROWS[i]] for i in S]))
            if best is not None and cost >= best[0]:
                continue
            if consistent(set(S)):
                best = (cost, [ROWS[i] for i in S])
        if best:
            break
    print(f'cheapest consistent set among price<=3 rows: {best}')
