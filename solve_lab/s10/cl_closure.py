"""CL step 5c: the FULL closure of the cluster system mod p, and the exact obstruction.

Iterate: columns = free inputs that move any row; rows = checks moved by any column.
Then Gaussian-eliminate [J | -r] over F_p with row-operation tracking, and verify any
inconsistent left-null vector y is a GLOBAL conserved functional (zero gradient over all
7,273 free inputs) via one reverse-AD pass.
"""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, fwdad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
BAD = [a for a in CHECKS if av0[a]]
FORB = {2081, 4287}

def grad_combo(weights):
    lam = collections.defaultdict(int)
    for a, wa in weights.items():
        if wa % P == 0: continue
        for w in L.avars[a]:
            lam[w] = (lam[w] + wa*ad.dpart(a, w, vm)) % P
    for t in reversed(ad.ORDER):
        lt = lam.get(t, 0) % P
        if not lt: continue
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0: continue
        f = -lt*pow(d, -1, P) % P
        for w in L.avars[a]:
            if w == t: continue
            dw = ad.dpart(a, w, vm)
            if dw: lam[w] = (lam[w] + f*dw) % P
        lam[t] = 0
    return {u: lam[u] % P for u in ad.FREESET if lam.get(u, 0) % P}

CACHE = os.path.join(HERE, 'cl_closure_cols.json')
if os.path.exists(CACHE):
    D = json.load(open(CACHE))
    cols = {int(u): {int(c): int(d) for c, d in mm.items()} for u, mm in D['cols'].items()}
    print(f'loaded {len(cols)} cached columns')
else:
    cols = {}
U = set()
for a in BAD: U |= set(ad.grad(a, vm))
U -= FORB
t0 = time.time()
for it in range(12):
    new = sorted(U - set(cols))
    for u in new:
        cols[u] = fwdad.jac_column(u, v0, vm, CHECKS)
    rows = sorted(set().union(*[set(c) for c in cols.values()]))
    # free inputs reachable from those rows
    U2 = set()
    for a in rows:
        U2 |= set(ad.grad(a, vm))
    U2 -= FORB
    print(f'  iter {it}: cols {len(cols)}  rows {len(rows)}  reachable free {len(U2)}  ({time.time()-t0:.0f}s)', flush=True)
    if U2 <= set(cols): break
    U = U2
json.dump({'cols': {str(u): {str(c): str(d) for c, d in m.items()} for u, m in cols.items()}},
          open(CACHE, 'w'))
U = sorted(cols)
rows = sorted(set().union(*[set(c) for c in cols.values()]) | set(BAD))
n, m = len(rows), len(U)
print(f'\nCLOSURE: {n} rows x {m} cols   ({time.time()-t0:.0f}s)')

ri = {a: i for i, a in enumerate(rows)}
A = [[0]*(m+1) for _ in range(n)]
for j, u in enumerate(U):
    for a, d in cols[u].items():
        A[ri[a]][j] = d % P
for a in rows:
    A[ri[a]][m] = (-av0[a]) % P
Tr = [{i: 1} for i in range(n)]

t0 = time.time()
r = 0
piv = []
for c in range(m):
    pr = None
    for i in range(r, n):
        if A[i][c]: pr = i; break
    if pr is None: continue
    A[r], A[pr] = A[pr], A[r]; Tr[r], Tr[pr] = Tr[pr], Tr[r]
    inv = pow(A[r][c], -1, P)
    if inv != 1:
        A[r] = [x*inv % P for x in A[r]]
        Tr[r] = {k: x*inv % P for k, x in Tr[r].items()}
    Ar, Trr = A[r], Tr[r]
    for i in range(r+1, n):
        f = A[i][c]
        if not f: continue
        Ai = A[i]
        A[i] = [(Ai[k] - f*Ar[k]) % P for k in range(m+1)]
        ti = Tr[i]
        for k, x in Trr.items():
            nvl = (ti.get(k, 0) - f*x) % P
            if nvl: ti[k] = nvl
            elif k in ti: del ti[k]
    piv.append(c); r += 1
    if r % 100 == 0: print(f'   pivot {r}/{m} ({time.time()-t0:.0f}s)', flush=True)
    if r == n: break
print(f'rank(J) = {r} of {m} columns, {n} rows   ({time.time()-t0:.0f}s)')

obs = [i for i in range(r, n) if A[i][m] % P and not any(A[i][j] for j in range(m))]
print(f'inconsistent left-null rows: {len(obs)}  (total left-null {n-r})')
if not obs:
    print('*** SYSTEM CONSISTENT: a mod-p linear repair of the whole cluster exists!')
    # back-substitute
    sol = [0]*m
    for i, c in enumerate(piv):
        s = A[i][m]
        for j in range(c+1, m):
            if A[i][j] and sol[j]: s = (s - A[i][j]*sol[j]) % P
        sol[c] = s % P
    json.dump({str(U[j]): str(sol[j]) for j in range(m) if sol[j]},
              open(os.path.join(HERE,'cl_delta.json'),'w'))
    print('saved cl_delta.json')
else:
    best = None
    for i in obs:
        y = {rows[k]: w for k, w in Tr[i].items() if w % P}
        g = grad_combo(y)
        val = sum(w*(av0[a] % P) for a, w in y.items()) % P
        print(f'  y: |support|={len(y)}  y.r={val}  global gradient support={len(g)}')
        if not g:
            print('  *** GLOBAL CONSERVED FUNCTIONAL -- certificate of infeasibility at this stratum')
            print(f'      support: {sorted(y)}')
            json.dump({'y': {str(a): str(w) for a, w in y.items()}, 'yr': str(val)},
                      open(os.path.join(HERE,'cl_cert.json'),'w'))
            best = y
            break
    if best is None:
        print('  no globally-conserved certificate among the inconsistent rows')
