"""S10 step 81: finish the budget-6 question, exactly.

The earlier exhaustive run timed out.  Reformulate so each test is trivial.

Dropping rows S leaves A_{-S} x = b_{-S}.  Its left null vectors, extended by
zeros on S, are exactly the y in leftnull(A) with supp(y) disjoint from S.  So

    consistent after dropping S  <=>  every y in leftnull(A) with y_S = 0 has y.b = 0
                                 <=>  t in colspace(Y[:, S])

where Y is a basis of leftnull(A) (49 x 128) and t = Y.b.  Each test is then a
49 x |S| rank check instead of a full elimination -- microseconds, not seconds.

This also yields the MINIMUM number of rows that must be sacrificed at all: the
fewest columns of Y whose span contains t.
"""
import os, sys, collections, itertools, time, json
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
FAIL = set(a for a in range(L.NA) if av[a] and a not in atom_out)
gcache = {}
def gradc(c):
    if c not in gcache:
        gcache[c] = {u: d for u, d in ad.grad(c, vm).items() if u not in BOOL}
    return gcache[c]
reach_cache = {}
def reach(u):
    if u in reach_cache: return reach_cache[u]
    seen, fr, ch = {u}, [u], set()
    while fr:
        nx = []
        for w in fr:
            for a in L.var_atoms[w]:
                if ad.dpart(a, w, vm) == 0: continue
                if a not in atom_out:
                    ch.add(a); continue
                t = atom_out[a][1]
                if t == w or t in seen or ad.dpart(a, t, vm) == 0: continue
                seen.add(t); nx.append(t)
        fr = nx
    reach_cache[u] = ch
    return ch

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
m = len(ROWS)
A = [[0] * n for _ in range(m)]
b = [0] * m
for i, c in enumerate(ROWS):
    for u, d in gradc(c).items():
        if u in cidx: A[i][cidx[u]] = d % P
    b[i] = (-av[c]) % P
print(f'closed system {m} rows x {n} cols')

# ---- left null space of A: rows y with y^T A = 0.  Work on A^T augmented. ----
# Solve A^T y = 0  (n equations, m unknowns)
M = [[A[i][j] for i in range(m)] for j in range(n)]      # n x m
piv = []; r = 0
for c in range(m):
    k = next((i for i in range(r, n) if M[i][c]), None)
    if k is None: continue
    M[r], M[k] = M[k], M[r]
    inv = pow(M[r][c], -1, P)
    M[r] = [x * inv % P for x in M[r]]
    for i in range(n):
        if i != r and M[i][c]:
            f = M[i][c]; Mr = M[r]
            M[i] = [(M[i][j] - f * Mr[j]) % P for j in range(m)]
    piv.append(c); r += 1
    if r == n: break
free = [c for c in range(m) if c not in piv]
Y = []
for fc in free:
    y = [0] * m
    y[fc] = 1
    for i, c in enumerate(piv):
        y[c] = (-M[i][fc]) % P
    Y.append(y)
print(f'rank(A) = {r}; left null space dimension = {len(Y)}')
t = [sum(Y[k][i] * b[i] for i in range(m)) % P for k in range(len(Y))]
print(f't = Y.b is {"ZERO (system consistent!)" if not any(t) else "nonzero -> inconsistent"}')

K = len(Y)
def in_span(S):
    """is t in the column span of Y[:, S] over GF(p)?"""
    cols_ = [[Y[k][i] for k in range(K)] for i in S]
    aug = [row[:] for row in cols_] + [t[:]]
    # row-reduce the |S| columns as rows of a |S| x K matrix, then reduce t
    basis = []
    for vec in cols_:
        w = vec[:]
        for bp, bv in basis:
            if w[bp]:
                f = w[bp] * pow(bv[bp], -1, P) % P
                w = [(w[j] - f * bv[j]) % P for j in range(K)]
        nz = next((j for j in range(K) if w[j]), None)
        if nz is not None:
            basis.append((nz, w))
    w = t[:]
    for bp, bv in basis:
        if w[bp]:
            f = w[bp] * pow(bv[bp], -1, P) % P
            w = [(w[j] - f * bv[j]) % P for j in range(K)]
    return not any(w)

price = {c: len(L.atom2eq.get(c, {})) for c in ROWS}
eqs = {c: set(L.atom2eq.get(c, ())) for c in ROWS}

t0 = time.time()
print('\n=== minimum number of rows that must be sacrificed (any price) ===')
minsize = None
for size in range(1, 6):
    hit = None
    for S in itertools.combinations(range(m), size):
        if in_span(S):
            hit = S; break
    print(f'  size {size}: {"FOUND " + str([ROWS[i] for i in hit]) if hit else "impossible"}'
          f'  ({time.time()-t0:.0f}s)', flush=True)
    if hit:
        minsize = (size, hit); break

print('\n=== exhaustive: any sacrifice set costing <= 6 equations? ===')
cheap = [i for i in range(m) if price[ROWS[i]] <= 6]
cheap.sort(key=lambda i: price[ROWS[i]])
print(f'rows priced <= 6: {len(cheap)}')
found = None
tested = [0]

def dfs(start, chosen, unioneq):
    """DFS with incremental cost pruning: abandon as soon as the union exceeds 6."""
    global found
    if found is not None:
        return
    if chosen:
        tested[0] += 1
        if in_span(tuple(chosen)):
            found = (len(unioneq), [ROWS[i] for i in chosen])
            return
    if len(chosen) == 6:
        return
    for k in range(start, len(cheap)):
        i = cheap[k]
        nu = unioneq | eqs[ROWS[i]]
        if len(nu) > 6:
            continue
        chosen.append(i)
        dfs(k + 1, chosen, nu)
        chosen.pop()
        if found is not None:
            return

dfs(0, [], set())
print(f'  cost-pruned DFS over all sizes 1..6: {tested[0]} within-budget sets tested, '
      f'{"FOUND" if found else "none"} ({time.time()-t0:.0f}s)', flush=True)
if found:
    print(f'\n*** WITHIN BUDGET: sacrifice {found[1]} at cost {found[0]} '
          f'=> score {L.NEQ - found[0]}')
    json.dump({'atoms': found[1], 'cost': found[0]},
              open(os.path.join(HERE, 'budget6fast.json'), 'w'))
else:
    print(f'\nEXHAUSTIVE NEGATIVE: no sacrifice set costs <= 6 equations '
          f'({tested} within-budget sets tested).')
    if minsize:
        sz, S = minsize
        print(f'  minimum sacrifice is {sz} rows: {[ROWS[i] for i in S]} '
              f'at cost {len(set().union(*[eqs[ROWS[i]] for i in S]))} equations')
