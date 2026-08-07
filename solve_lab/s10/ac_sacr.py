"""S12 step 7: the closure as a MINIMUM-EQUATION-COST decoding problem.

The closure M (rows = checks, cols = free inputs) has full column rank, so the
solve is *injective*: choosing which rows to pin fixes x, and the violated set is
D = supp(Mx - r).  The final score is 39033 - |eqs_of_atoms(D)|.  x = 0 gives the
current D = {21617,29539,37662,40826} costing 24 equations.  Minimise over the
choice of pinned rows.
"""
import os, sys, json, time, random, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FORBID = {2081, 4287}
CHECKS = sorted(a for a in range(L.NA) if a not in atom_out)

def closure_of(v, BAD, iters=9):
    vm = [x % P for x in v]
    gc = {}
    def gr(c):
        if c not in gc: gc[c] = set(ad.grad(c, vm)) - FORBID
        return gc[c]
    U, rows, cols = set(), set(BAD), {}
    for it in range(iters):
        newU = set()
        for c in rows: newU |= gr(c)
        newU -= U; U |= newU
        for u in sorted(newU): cols[u] = jac_column(u, v, vm, CHECKS)
        nr = set(BAD)
        for u in U: nr |= set(cols[u])
        grew = nr - rows; rows = nr
        if not newU and not grew: break
    return sorted(rows), sorted(U), cols

def decode(rows, Us, cols, av, order):
    """Gauss-Jordan with rows presented in `order`; returns x and violated set."""
    perm = order
    ri = {c: i for i, c in enumerate(perm)}
    n, m = len(perm), len(Us)
    M = [[0]*(m+1) for _ in range(n)]
    for j, u in enumerate(Us):
        for c, d in cols[u].items():
            if c in ri: M[ri[c]][j] = d % P
    for c in perm: M[ri[c]][m] = (-av[c]) % P
    r_ = 0; piv = []
    for j in range(m):
        k = next((i for i in range(r_, n) if M[i][j]), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]
        inv = pow(M[r_][j], -1, P)
        M[r_] = [x*inv % P for x in M[r_]]
        for i in range(n):
            if i != r_ and M[i][j]:
                f = M[i][j]
                M[i] = [(a2-f*b2) % P for a2, b2 in zip(M[i], M[r_])]
        piv.append(j); r_ += 1
    x = [0]*m
    for k, j in enumerate(piv): x[j] = M[k][m]
    nbad = sum(1 for i in range(r_, n) if M[i][m])
    # residual rho = M x - r, evaluated sparsely
    rho = collections.defaultdict(int)
    for c in perm: rho[c] = av[c] % P
    for j, u in enumerate(Us):
        if not x[j]: continue
        for c, d in cols[u].items():
            if c in ri: rho[c] = (rho[c] + d*x[j]) % P
    D = sorted(c for c in perm if rho[c] % P)
    return x, D, r_, nbad

if __name__ == '__main__':
    v = L.load(os.path.join(HERE,'mod9118_0.json'))
    av = L.all_atom_values(v)
    BAD = [21617, 29539]
    t0 = time.time()
    rows, Us, cols = closure_of(v, BAD)
    print(f'closure {len(rows)} x {len(Us)} ({time.time()-t0:.0f}s)', flush=True)
    w = {c: len(L.atom2eq.get(c, {})) for c in rows}
    print(f'row weights (equations per check): '
          f'{dict(sorted(collections.Counter(w.values()).items()))}')
    base = sorted(L.eqs_of_atoms([a for a in range(L.NA) if av[a]]))
    print(f'baseline D = nonzero atoms, cost {len(base)} equations -> score {L.NEQ-len(base)}')
    random.seed(5)
    results = []
    orders = {'natural': rows,
              'cost-desc': sorted(rows, key=lambda c: (-w[c], c)),
              'cost-asc': sorted(rows, key=lambda c: (w[c], c))}
    for k in range(6):
        o = list(rows); random.shuffle(o); orders[f'rand{k}'] = o
    for name, o in orders.items():
        t1 = time.time()
        x, D, rank, nbad = decode(rows, Us, cols, av, o)
        eqs = sorted(L.eqs_of_atoms(D))
        results.append((len(eqs), name, D))
        print(f'  {name:>10}: rank {rank}/{len(Us)}  violated rows {len(D)}  '
              f'equation cost {len(eqs)}  -> score {L.NEQ-len(eqs)}  '
              f'({time.time()-t1:.0f}s)', flush=True)
        print(f'             D = {D}')
    results.sort()
    json.dump({'best_cost': results[0][0], 'best_D': results[0][2],
               'rows': rows, 'Us': Us}, open(os.path.join(HERE,'ac_sacr.json'),'w'))
    print(f'\nBEST equation cost over the orderings tried: {results[0][0]} '
          f'({results[0][1]}) -> mod-p ceiling score {L.NEQ-results[0][0]}')
