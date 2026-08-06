"""S10 step 34: Newton over GF(p) on the check system.

Build J = d(failing checks)/d(free inputs) mod p by reverse-mode AD, solve
J.delta = -residual over GF(p) restricted to NON-boolean free inputs, apply,
forward-evaluate, repeat.  Prefer free inputs with small footprints so the step
disturbs as little as possible.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out

# ---- identify boolean free inputs (some atom pins x*x - x) -----------------
BOOL = set()
for a, poly in enumerate(L.polys):
    ks = list(poly.items())
    if len(ks) == 2:
        sq = [m for m, c in ks if len(m) == 2 and m[0] == m[1]]
        li = [m for m, c in ks if len(m) == 1]
        if sq and li and sq[0][0] == li[0][0]:
            BOOL.add(li[0][0])
print(f'boolean variables detected: {len(BOOL)}')
CAND = [u for u in ad.FREE if u not in BOOL]
print(f'non-boolean free inputs: {len(CAND)}')


def solve_gfp(rows, rhs, cols):
    """Any solution x (dict col->val) of rows.x = rhs over GF(p), or None."""
    m = len(rows)
    A = [[r.get(c, 0) % P for c in cols] + [rhs[i] % P] for i, r in enumerate(rows)]
    n = len(cols)
    piv = []
    r = 0
    for c in range(n):
        k = next((i for i in range(r, m) if A[i][c] % P), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        inv = pow(A[r][c], -1, P)
        A[r] = [x * inv % P for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % P for j in range(n + 1)]
        piv.append(c); r += 1
        if r == m:
            break
    for i in range(r, m):
        if A[i][n] % P:
            return None                      # inconsistent
    x = {}
    for i, c in enumerate(piv):
        x[cols[c]] = A[i][n] % P
    return x


def status(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    return av, nz, fail


v = L.load(os.path.join(HERE, 'forward_state.json'))
best = None
for it in range(25):
    av, nz, fail = status(v)
    sc = L.NEQ - len(fail)
    checks = [a for a in nz if a not in atom_out]
    gates = [a for a in nz if a in atom_out]
    print(f'\niter {it}: nz={len(nz)} (checks {len(checks)}, gates {len(gates)}) '
          f'failing={len(fail)} score={sc}', flush=True)
    print(f'   checks: {checks}', flush=True)
    if best is None or sc > best[0]:
        best = (sc, list(v))
    if not nz:
        print('*** ALL ATOMS ZERO -- FULL SOLUTION ***', flush=True)
        break
    if not checks:
        print('   only gate atoms left', flush=True); break
    vm = [x % P for x in v]
    grads = []
    support = set()
    for c in checks:
        g = ad.grad(c, vm)
        g = {u: d for u, d in g.items() if u not in BOOL}
        grads.append(g); support |= set(g)
    cols = sorted(support, key=lambda u: (len(L.var_atoms[u]), len(L.var_eqs[u])))
    print(f'   support {len(cols)} non-boolean free inputs', flush=True)
    rhs = [(-av[c]) % P for c in checks]
    sol = solve_gfp(grads, rhs, cols)
    if sol is None:
        print('   *** GF(p) system INCONSISTENT at this point', flush=True)
        break
    print(f'   step moves {len(sol)} inputs: '
          f'{[f"x_{u}" for u in list(sol)[:8]]}', flush=True)
    for u, d in sol.items():
        v[u] = v[u] + d
    ad.fwd(v)

av, nz, fail = status(v)
print(f'\nFINAL nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}')
print(f'BEST {best[0]}')
T.save(best[1], os.path.join(HERE, 'newton_best.json'))
