"""S10 step 58: extract the inconsistency CERTIFICATES, then min-cost hitting set.

Augment [A | b | I] and eliminate on A's columns.  Any row whose A-part is zero
but whose b-part is not gives a left-null vector y with y.A = 0, y.b != 0 -- an
explicit certificate of inconsistency, and its I-part names exactly which original
checks combine to produce it.

To make the system solvable we must SACRIFICE a set of checks hitting every
certificate.  Minimise the equation cost of that set; the budget to beat is 7.
"""
import os, sys, collections, json, time
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
n, m = len(COLS), len(ROWS)
cidx = {u: i for i, u in enumerate(COLS)}
print(f'closed system {m} rows x {n} cols', flush=True)

# [A | b | I]
M = []
for i, c in enumerate(ROWS):
    r = [0] * (n + 1 + m)
    for u, d in gradc(c).items():
        if u in cidx: r[cidx[u]] = d % P
    r[n] = (-av[c]) % P
    r[n + 1 + i] = 1
    M.append(r)

rr = 0
for col in range(n):
    k = next((i for i in range(rr, m) if M[i][col]), None)
    if k is None: continue
    M[rr], M[k] = M[k], M[rr]
    inv = pow(M[rr][col], -1, P)
    M[rr] = [x * inv % P for x in M[rr]]
    for i in range(m):
        if i != rr and M[i][col]:
            f = M[i][col]; Mr = M[rr]
            M[i] = [(M[i][j] - f * Mr[j]) % P for j in range(n + 1 + m)]
    rr += 1
    if rr == m: break
print(f'rank(A) = {rr}', flush=True)

certs = []
for i in range(rr, m):
    if M[i][n]:
        sup = [ROWS[j] for j in range(m) if M[i][n + 1 + j]]
        certs.append(sup)
print(f'\nINCONSISTENCY CERTIFICATES: {len(certs)}', flush=True)
price = {c: len(L.atom2eq.get(c, {})) for c in ROWS}
for k, s in enumerate(certs):
    cheapest = sorted(((price[c], c) for c in s))[:6]
    print(f'  cert {k}: support {len(s)} checks; cheapest members {cheapest}', flush=True)

# min-cost hitting set (greedy + swap improvement)
def cost_of(S):
    eqs = set()
    for c in S: eqs |= set(L.atom2eq.get(c, ()))
    return len(eqs)

uncovered = [set(s) for s in certs]
S = set()
while uncovered:
    bestc, bestv = None, None
    for c in set().union(*uncovered):
        hits = sum(1 for u in uncovered if c in u)
        val = (cost_of(S | {c}) - cost_of(S)) / hits if hits else 1e9
        if bestv is None or val < bestv:
            bestv, bestc = val, c
    S.add(bestc)
    uncovered = [u for u in uncovered if bestc not in u]
print(f'\ngreedy hitting set: {sorted(S)}')
print(f'  prices: {[(c, price[c]) for c in sorted(S)]}')
print(f'  TOTAL EQUATION COST = {cost_of(S)}   (budget to beat: 7)')
json.dump({'certs': certs, 'S': sorted(S), 'cost': cost_of(S),
           'price': {str(k): v for k, v in price.items()}},
          open(os.path.join(HERE, 'certs.json'), 'w'))
