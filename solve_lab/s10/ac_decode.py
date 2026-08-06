"""S12 step 9: proper minimum-EQUATION-COST decoding of the closure coset.

The score after a linearised move is 39033 - |eqs_of_atoms(D)|, D = supp(Mx - r).
396 of the 1655 closure checks touch only ONE equation.  So ask directly: is the
subsystem on the EXPENSIVE rows consistent?  If it is, the whole cost can be
pushed onto weight-1 checks and the ceiling is the number of distinct equations
they touch.
"""
import os, sys, json, time, collections, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from ac_sacr import closure_of
P = ad.P

def build(rows, Us, cols, av):
    ri = {c: i for i, c in enumerate(rows)}
    n, m = len(rows), len(Us)
    M = [[0]*(m+1) for _ in range(n)]
    for j, u in enumerate(Us):
        for c, d in cols[u].items():
            if c in ri: M[ri[c]][j] = d % P
    for c in rows: M[ri[c]][m] = (-av[c]) % P
    return M, n, m

def gj(M, n, m):
    """Gauss-Jordan on the augmented matrix; returns (rank, piv, nbad, M)."""
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
                M[i] = [(a-f*b) % P for a, b in zip(M[i], M[r_])]
        piv.append(j); r_ += 1
    nbad = sum(1 for i in range(r_, n) if M[i][m])
    return r_, piv, nbad, M

def resid(rows, Us, cols, av, x):
    rho = {c: av[c] % P for c in rows}
    ri = set(rows)
    for j, u in enumerate(Us):
        if not x[j]: continue
        for c, d in cols[u].items():
            if c in ri: rho[c] = (rho[c] + d*x[j]) % P
    return sorted(c for c in rows if rho[c] % P)

if __name__ == '__main__':
    v = L.load(os.path.join(HERE,'mod9118_0.json'))
    av = L.all_atom_values(v)
    BAD = [21617, 29539]
    t0 = time.time()
    rows, Us, cols = closure_of(v, BAD)
    print(f'closure {len(rows)} x {len(Us)}  ({time.time()-t0:.0f}s)', flush=True)
    w = {c: len(L.atom2eq.get(c, {})) for c in rows}
    for thr in (2, 3, 5, 8, 10, 12):
        sub = [c for c in rows if w[c] >= thr]
        M, n, m = build(sub, Us, cols, av)
        rank, piv, nbad, M = gj(M, n, m)
        x = [0]*m
        for k, j in enumerate(piv): x[j] = M[k][m]
        ok = (nbad == 0)
        note = ''
        if ok:
            D = resid(rows, Us, cols, av, x)
            eqs = L.eqs_of_atoms(D)
            note = (f'  CONSISTENT -> violated {len(D)} cheap checks, '
                    f'{len(eqs)} equations, score {L.NEQ-len(eqs)}')
        print(f'  rows with weight >= {thr:>2}: {n:>5} rows, rank {rank}/{m}, '
              f'kernel {m-rank}, inconsistent {nbad}{note}', flush=True)
    # also: which cheap rows can absorb the inconsistency?
    print(f'\nweight-1 checks in the closure: {sum(1 for c in rows if w[c]==1)}; '
          f'they touch {len(L.eqs_of_atoms([c for c in rows if w[c]==1]))} distinct equations')
