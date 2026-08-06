"""CL step 5b: solve J*delta = -r over F_p; extract the left-kernel obstruction and
verify it is a GLOBAL conserved functional by one reverse-AD pass."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm = [x % P for x in v0]
av0 = L.all_atom_values(v0)

# ---------- generalised reverse AD: gradient of sum_a w_a * atom_a ----------
def grad_combo(weights, vm):
    lam = collections.defaultdict(int)
    for a, wa in weights.items():
        if wa % P == 0: continue
        for w in L.avars[a]:
            lam[w] = (lam[w] + wa * ad.dpart(a, w, vm)) % P
    for t in reversed(ad.ORDER):
        lt = lam.get(t, 0) % P
        if not lt: continue
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0: continue
        f = -lt * pow(d, -1, P) % P
        for w in L.avars[a]:
            if w == t: continue
            dw = ad.dpart(a, w, vm)
            if dw: lam[w] = (lam[w] + f * dw) % P
        lam[t] = 0
    return {u: lam[u] % P for u in ad.FREESET if lam.get(u, 0) % P}

D = json.load(open(os.path.join(HERE,'cl_jac.json')))
U = D['U']; rows = D['rows']
cols = {int(u): {int(c): int(d) for c, d in m.items()} for u, m in D['cols'].items()}
resid = {int(r): int(x) for r, x in D['resid'].items()}
ri = {a: i for i, a in enumerate(rows)}
n, m = len(rows), len(U)
print(f'system {n} rows x {m} cols')

# augmented matrix M = [J | -r] with an identity block to track row ops -> left kernel
A = [[0]*(m+1) for _ in range(n)]
for j, u in enumerate(U):
    for a, d in cols[u].items():
        A[ri[a]][j] = d % P
for a, x in resid.items():
    A[ri[a]][m] = (-x) % P
Tr = [[1 if i == k else 0 for k in range(n)] for i in range(n)]   # row-op record

def elim(ncols):
    piv = []
    r = 0
    for c in range(ncols):
        pr = None
        for i in range(r, n):
            if A[i][c] % P: pr = i; break
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; Tr[r], Tr[pr] = Tr[pr], Tr[r]
        inv = pow(A[r][c], -1, P)
        A[r] = [x*inv % P for x in A[r]]; Tr[r] = [x*inv % P for x in Tr[r]]
        for i in range(n):
            if i != r and A[i][c] % P:
                f = A[i][c]
                A[i] = [(A[i][k] - f*A[r][k]) % P for k in range(m+1)]
                Tr[i] = [(Tr[i][k] - f*Tr[r][k]) % P for k in range(n)]
        piv.append(c); r += 1
        if r == n: break
    return piv, r

t0 = time.time()
piv, rk = elim(m)
print(f'rank(J) = {rk}  (of {m} columns)   [{time.time()-t0:.0f}s]')
# rows r..n-1 are now all-zero in the J block: their last entry must be 0 for consistency
obstr = []
for i in range(rk, n):
    if any(A[i][j] % P for j in range(m)):
        continue
    if A[i][m] % P:
        obstr.append(i)
print(f'inconsistent left-null rows: {len(obstr)} (of {n-rk} left-null vectors)')

if not obstr:
    print('SYSTEM IS CONSISTENT -- a linear repair exists!')
    sol = [0]*m
    for i, c in enumerate(piv):
        sol[c] = A[i][m] % P
    json.dump({'U': U, 'delta': {str(U[j]): str(sol[j]) for j in range(m) if sol[j]}},
              open(os.path.join(HERE,'cl_delta.json'),'w'))
    print('saved cl_delta.json')
else:
    print('\n=== obstruction functionals ===')
    for i in obstr[:5]:
        y = {rows[k]: Tr[i][k] % P for k in range(n) if Tr[i][k] % P}
        val = sum(w * (av0[a] % P) for a, w in y.items()) % P
        print(f'  y supported on {len(y)} checks; y.r = {val}')
        top = sorted(y.items(), key=lambda kv: -abs(av0[kv[0]] % P))[:12]
        print(f'    checks with nonzero residue in support: '
              f'{[(f"a{a}", str(w)[:14]) for a, w in y.items() if av0[a]][:8]}')
        # GLOBAL verification: is grad(sum y_a a_a) identically zero over ALL free inputs?
        g = grad_combo(y, vm)
        print(f'    GLOBAL gradient support of the functional: {len(g)} free inputs '
              f'{"-> CONSERVED (certificate)" if not g else "-> not conserved; escapes:"}')
        if g:
            print(f'    escaping free inputs (first 20): {sorted(g)[:20]}')
        json.dump({'y': {str(a): str(w) for a, w in y.items()},
                   'yr': str(val), 'escapes': sorted(g)},
                  open(os.path.join(HERE, f'cl_obstr_{i}.json'), 'w'))
