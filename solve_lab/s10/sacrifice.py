"""S10 step 56: BUDGET ATTACK.

The current branch pays 7 failing equations, so any structural violation costing
<= 6 equations beats it.  The closed constrained system has exactly 6 independent
inconsistencies, so the question is not "can everything be satisfied" (it cannot)
but:

    which minimal-cost set of checks do I SACRIFICE so the rest becomes solvable?

cost(check) = number of equations that check lives in.  Square checks and the two
big checks cost 1 each -- the cheapest guards in the instance.
"""
import os, sys, collections, json, time, itertools
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
print('failing checks:', sorted(FAIL), flush=True)

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
    rows_set = set()
    for u in cols: rows_set |= reach(u)
    new = set(cols)
    for c in rows_set: new |= set(gradc(c))
    if len(new) == len(cols): break
    cols = new
rows_set = set()
for u in cols: rows_set |= reach(u)
COLS = sorted(cols)
ROWS = []
for c in sorted(rows_set):
    g = {u: d for u, d in gradc(c).items() if u in cols}
    if not g and (av[c] % P) == 0: continue
    ROWS.append(c)
print(f'closed system: {len(ROWS)} rows x {len(COLS)} cols', flush=True)

cost = {c: len(L.atom2eq.get(c, {})) for c in ROWS}
print('cheapest rows by equation cost:',
      sorted(((cost[c], c) for c in ROWS))[:18], flush=True)


def consistent(drop):
    keep = [c for c in ROWS if c not in drop]
    n = len(COLS); cidx = {u: i for i, u in enumerate(COLS)}
    M = []
    for c in keep:
        r = [0] * (n + 1)
        for u, d in gradc(c).items():
            if u in cidx: r[cidx[u]] = d % P
        r[n] = (-av[c]) % P
        M.append(r)
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


t0 = time.time()
print(f'\nfull system consistent? {consistent(set())}', flush=True)

# greedy: repeatedly drop the cheapest row that most reduces inconsistency
cand = sorted(ROWS, key=lambda c: cost[c])
best = None
for k in range(1, 7):
    found = False
    for combo in itertools.combinations(cand[:26], k):
        if sum(cost[c] for c in combo) > 6:
            continue
        if consistent(set(combo)):
            print(f'\n*** CONSISTENT after sacrificing {combo} '
                  f'(equation cost {sum(cost[c] for c in combo)})', flush=True)
            best = combo; found = True
            break
    print(f'  size {k}: {"FOUND" if found else "none within budget"} '
          f'({time.time()-t0:.0f}s)', flush=True)
    if found:
        break
if not best:
    print('\nno sacrifice set within the 6-equation budget among the cheapest rows')
    # what is the cheapest consistent sacrifice at all?
    for k in range(1, 5):
        got = None
        for combo in itertools.combinations(cand[:20], k):
            if consistent(set(combo)):
                got = combo; break
        if got:
            print(f'  cheapest consistent sacrifice of size {k}: {got} '
                  f'cost={sum(cost[c] for c in got)}')
            break
