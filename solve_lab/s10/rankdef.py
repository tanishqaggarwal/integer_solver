"""S10 step 39: how BIG is the obstruction?

One Gaussian elimination on [A | b] over GF(p) gives:
  rank(A), rank([A|b]), and the number of independent inconsistencies
  = rank([A|b]) - rank(A).

If that number is small, dropping that many checks makes the rest solvable, and
each dropped check costs a known number of equations -- which converts the
obstruction directly into a score bound.
"""
import os, sys, json, time, collections
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
USE_BOOL = os.environ.get('NOBOOL', '') == ''
STATE = os.environ.get('STATE', os.path.join(HERE, 'forward_state.json'))

v = L.load(STATE)
vm = [x % P for x in v]
av = L.all_atom_values(v)
FAIL = set(a for a in range(L.NA) if av[a] and a not in atom_out)
print(f'state={os.path.basename(STATE)} failing checks={sorted(FAIL)} bits={"in" if USE_BOOL else "out"}',
      flush=True)

gcache = {}
def gradc(c):
    if c not in gcache:
        g = ad.grad(c, vm)
        if not USE_BOOL:
            g = {u: d for u, d in g.items() if u not in BOOL}
        gcache[c] = g
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
for rnd in range(10):
    rows_set = set()
    for u in cols: rows_set |= reach(u)
    newcols = set(cols)
    for c in rows_set: newcols |= set(gradc(c))
    if len(newcols) == len(cols): break
    cols = newcols
rows_set = set()
for u in cols: rows_set |= reach(u)
COLS = sorted(cols); ROWS = sorted(rows_set)
print(f'closed system {len(ROWS)} x {len(COLS)}', flush=True)

R, RHS, NAME = [], [], []
for c in ROWS:
    g = {u: d for u, d in gradc(c).items() if u in cols}
    if not g and (av[c] % P) == 0: continue
    R.append(g); RHS.append((-av[c]) % P); NAME.append(c)
n = len(COLS); cidx = {c: i for i, c in enumerate(COLS)}
M = [[0] * (n + 1) for _ in R]
for i, g in enumerate(R):
    for u, d in g.items(): M[i][cidx[u]] = d % P
    M[i][n] = RHS[i] % P
print(f'matrix {len(M)} x {n+1}', flush=True)

t0 = time.time()
r = 0; piv = []
bad = []
for c in range(n):
    k = next((i for i in range(r, len(M)) if M[i][c]), None)
    if k is None: continue
    M[r], M[k] = M[k], M[r]
    NAME[r], NAME[k] = NAME[k], NAME[r]
    inv = pow(M[r][c], -1, P)
    M[r] = [x * inv % P for x in M[r]]
    for i in range(len(M)):
        if i != r and M[i][c]:
            f = M[i][c]
            Mr = M[r]
            M[i] = [(M[i][j] - f * Mr[j]) % P for j in range(n + 1)]
    piv.append(c); r += 1
    if r == len(M): break
rankA = r
inc = [i for i in range(r, len(M)) if M[i][n]]
print(f'\nrank(A)      = {rankA}')
print(f'rows         = {len(M)}, cols = {n}')
print(f'INDEPENDENT INCONSISTENCIES = {len(inc)}   ({time.time()-t0:.0f}s)')
print(f'  -> at least {len(inc)} check(s) cannot be satisfied by any first-order move')
if inc:
    print(f'  witness rows (post-elimination indices): {inc[:10]}')
