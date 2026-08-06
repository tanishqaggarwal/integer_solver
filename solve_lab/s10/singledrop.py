"""S10 step 75: the last door -- what is the MINIMUM-COST sacrifice, exactly?

Correction to Part III: b is a single column, so rank([A|b]) - rank(A) <= 1.
There is ONE independent inconsistency, not six; the "6 certificates" were six
witnessing rows of the same obstruction.  So a SINGLE dropped row may suffice.

Test every one of the 128 rows: drop it, re-check consistency over GF(p), and
report the price of every row that works.  Then all pairs among the cheap rows.
Budget to beat: 7.
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
print(f'closed system {len(ROWS)} rows x {n} cols', flush=True)

FULL = []
for c in ROWS:
    r = [0] * (n + 1)
    for u, d in gradc(c).items():
        if u in cidx: r[cidx[u]] = d % P
    r[n] = (-av[c]) % P
    FULL.append(r)


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


price = {c: len(L.atom2eq.get(c, {})) for c in ROWS}
print(f'full system consistent: {consistent(set())}', flush=True)

t0 = time.time()
print('\n=== every single-row drop ===', flush=True)
works = []
for i, c in enumerate(ROWS):
    if consistent({i}):
        works.append((price[c], c, i))
        print(f'  DROP a{c} (price {price[c]}) -> CONSISTENT', flush=True)
works.sort()
print(f'single drops that work: {len(works)}  ({time.time()-t0:.0f}s)')
if works:
    print(f'  cheapest: a{works[0][1]} at price {works[0][0]}  '
          f'=> score {L.NEQ - works[0][0]}')
    json.dump({'works': [[p, c] for p, c, i in works]},
              open(os.path.join(HERE, 'singledrop.json'), 'w'))
else:
    print('  none. trying pairs among the 30 cheapest rows...', flush=True)
    cheap = sorted(range(len(ROWS)), key=lambda i: price[ROWS[i]])[:30]
    best = None
    for i, j in itertools.combinations(cheap, 2):
        cost = len(set(L.atom2eq.get(ROWS[i], ())) | set(L.atom2eq.get(ROWS[j], ())))
        if best is not None and cost >= best[0]:
            continue
        if consistent({i, j}):
            best = (cost, ROWS[i], ROWS[j])
            print(f'  PAIR a{ROWS[i]}+a{ROWS[j]} cost {cost} -> CONSISTENT', flush=True)
    print(f'best pair: {best}  ({time.time()-t0:.0f}s)')
