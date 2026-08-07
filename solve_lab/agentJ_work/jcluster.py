#!/usr/bin/env python3
"""The deliverable's own cluster: 12 equations, 7 atom values.

score = 39033 - (12 - k)  where k = number of the 12 rows driven to zero.
Deliverable => k = 5.  k = 6 would be 39027.  k = 12 would be a full solve.

Reachable set: start at the deliverable x*, find every variable whose perturbation
keeps the nonzero-atom set inside T (so nothing else breaks), and record the exact
induced change in the 7 atom values.  That generates a lattice L; the reachable
atom vectors are a* + L.  Then for every subset S of the 12 rows ask whether
C_S (a* + L t) = 0 has an INTEGER solution t (Hermite normal form).
"""
import os, pickle, sys, itertools, json
from collections import defaultdict
import jengine as E

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
eqs, atoms, polys = M['eqs'], M['atoms'], P['polys']
NV = E.NV
occ = defaultdict(list)
for i, s in enumerate(E.varsof):
    for v in s:
        occ[v].append(i)

DEL = os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json')


def atomval(i, val):
    s = 0
    for k, c in polys[i].items():
        t = c
        for j in k:
            t *= val[j]
        s += t
    return s


def all_nonzero(val):
    return [i for i in range(len(polys)) if atomval(i, val) != 0]


# ---------- integer linear algebra: solve B t = c over Z ----------
def hnf_solve(B, c):
    """B: list of rows (lists) m x n ; c: list length m.  Return integer t or None."""
    m = len(B); n = len(B[0]) if m else 0
    A = [row[:] for row in B]
    rhs = c[:]
    V = [[1 if i == j else 0 for j in range(n)] for i in range(n)]  # column ops
    row = 0
    piv = []
    for col in range(n):
        if row >= m:
            break
        # find a row >= `row` with nonzero entry in some column >= col; use col-gcd
        # gcd elimination across columns col..n-1 within row `row`
        # first make sure row `row` has a nonzero entry in columns >= col
        r = row
        while r < m and all(A[r][cc] == 0 for cc in range(col, n)):
            r += 1
        if r == m:
            break
        if r != row:
            A[r], A[row] = A[row], A[r]
            rhs[r], rhs[row] = rhs[row], rhs[r]
        # column-reduce row `row` to a single nonzero at `col`
        while True:
            nz = [cc for cc in range(col, n) if A[row][cc] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda cc: abs(A[row][cc]))
            c0 = nz[0]
            for cc in nz[1:]:
                q = A[row][cc] // A[row][c0]
                if q:
                    for rr in range(m):
                        A[rr][cc] -= q * A[rr][c0]
                    for rr in range(n):
                        V[rr][cc] -= q * V[rr][c0]
            # re-loop
        nz = [cc for cc in range(col, n) if A[row][cc] != 0]
        if not nz:
            row += 1
            continue
        c0 = nz[0]
        if c0 != col:
            for rr in range(m):
                A[rr][col], A[rr][c0] = A[rr][c0], A[rr][col]
            for rr in range(n):
                V[rr][col], V[rr][c0] = V[rr][c0], V[rr][col]
        piv.append((row, col))
        row += 1
    # now A is lower-trapezoidal in the pivot structure: solve by forward substitution
    y = [0] * n
    used = set()
    for (r, cc) in piv:
        s = rhs[r]
        for k in range(n):
            if k != cc and y[k]:
                s -= A[r][k] * y[k]
        if A[r][cc] == 0:
            if s != 0:
                return None
            continue
        if s % A[r][cc] != 0:
            return None
        y[cc] = s // A[r][cc]
        used.add(cc)
    # verify all rows
    for r in range(m):
        s = sum(A[r][k] * y[k] for k in range(n))
        if s != rhs[r]:
            return None
    t = [sum(V[i][k] * y[k] for k in range(n)) for i in range(n)]
    return t


if __name__ == '__main__':
    val = E.load(DEL)
    T = all_nonzero(val)
    print("deliverable nonzero atoms T =", T)
    a0 = [atomval(i, val) for i in T]

    # the equations touching T, and the coefficient matrix
    rows = []
    reqs = []
    for e in eqs:
        row = {}
        for c, j in e['terms']:
            row[j] = row.get(j, 0) + c
        if any(j in T and row[j] for j in row):
            reqs.append(e['i'])
            rows.append([row.get(j, 0) for j in T])
    print(f"equations touching T: {len(reqs)} -> {reqs}")
    cur = [sum(rows[r][k] * a0[k] for k in range(len(T))) for r in range(len(rows))]
    print("currently zero rows:", sum(1 for x in cur if x == 0), "of", len(rows))

    # ---- generators: variables all of whose atoms lie in T ----
    Tset = set(T)
    gens = []
    gnames = []
    varsT = set()
    for i in T:
        varsT |= E.varsof[i]
    for v in sorted(varsT):
        if set(occ[v]) <= Tset:
            v2 = list(val); v2[v] += 1
            d = [atomval(i, v2) - atomval(i, val) for i in T]
            if any(d):
                gens.append(d); gnames.append(f"x_{v}")
    print(f"\nstrictly-safe generators ({len(gens)}): {gnames}")
    for n_, g in zip(gnames, gens):
        print("   ", n_, g)

    # ---- the subset search ----
    L = gens
    nT = len(T)
    print("\nsubset search: can k of the %d rows be driven to zero?" % len(rows))
    best = None
    for k in range(len(rows), 4, -1):
        found = None
        for S in itertools.combinations(range(len(rows)), k):
            B = [[sum(rows[r][c] * L[g][c] for c in range(nT)) for g in range(len(L))]
                 for r in S]
            c = [-sum(rows[r][c] * a0[c] for c in range(nT)) for r in S]
            t = hnf_solve(B, c)
            if t is not None:
                found = (S, t)
                break
        if found:
            print(f"  k = {k}: SOLVABLE  rows {found[0]}  => score {39033 - (len(rows)-k)}")
            best = (k, found)
            break
        else:
            print(f"  k = {k}: no integral solution over {len(list(itertools.combinations(range(len(rows)),k)))} subsets")
    pickle.dump({'T': T, 'reqs': reqs, 'rows': rows, 'a0': a0, 'gens': gens,
                 'gnames': gnames, 'best': best},
                open(os.path.join(HERE, 'jcluster.pkl'), 'wb'))
