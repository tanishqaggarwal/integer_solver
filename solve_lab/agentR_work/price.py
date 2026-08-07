#!/usr/bin/env python3
"""Price a configuration's defect footprint the way the lab priced the deliverable's:
build the (equations touched) x (nonzero atoms) integer coefficient matrix and ask how many
of those equations a nonzero atom vector can kill.  Upper bound on the achievable score."""
import sys, os, json, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import run_cfg, E
from fractions import Fraction

def price(bits):
    sc, nz, ok, v = run_cfg(bits)
    r = E.run(v)
    nzi = [i for i, x in enumerate(r) if x]
    S = set(nzi)
    rows = []
    for j, rr in enumerate(E.eqres):
        c = {i: k for k, i in rr if i in S}
        if c: rows.append((j, [c.get(i, 0) for i in nzi]))
    M = [row for _, row in rows]
    m, n = len(M), len(nzi)
    # max rows killable by a nonzero vector a: try null vectors of every (n-1)-subset of rows
    best, bestvec = 0, None
    import random
    from itertools import combinations
    def nullvec(sub):
        # exact rational nullspace of the sub matrix (n columns)
        A = [list(map(Fraction, x)) for x in sub]
        piv = []; rr = 0
        for c in range(n):
            p = next((i for i in range(rr, len(A)) if A[i][c] != 0), None)
            if p is None: continue
            A[rr], A[p] = A[p], A[rr]
            pv = A[rr][c]
            A[rr] = [x / pv for x in A[rr]]
            for i in range(len(A)):
                if i != rr and A[i][c] != 0:
                    f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[rr])]
            piv.append(c); rr += 1
            if rr == len(A): break
        free = [c for c in range(n) if c not in piv]
        if not free: return None
        f0 = free[0]
        x = [Fraction(0)] * n; x[f0] = Fraction(1)
        for i, c in enumerate(piv): x[c] = -A[i][f0]
        return x
    cand = []
    if n >= 2:
        for sub in combinations(range(m), min(n - 1, m)):
            x = nullvec([M[i] for i in sub])
            if x: cand.append(x)
    seen = set()
    for x in cand:
        key = tuple(x)
        if key in seen: continue
        seen.add(key)
        cnt = sum(1 for row in M if sum(a * b for a, b in zip(row, x)) == 0)
        if cnt > best: best, bestvec = cnt, x
    return dict(bits=bits, score=sc, atoms=n, eqs=m, max_killable=best,
                floor_failing=m - best, implied_score=39033 - (m - best))

if __name__ == '__main__':
    out = {}
    for c in ([24601, 2081], [24601], [47], [91], [542]):
        p = price(c)
        out[','.join(map(str, c))] = p
        print('%-16s atoms=%d eqs=%d  max killable=%d -> floor failing=%d -> implied score <= %d'
              % (str(c), p['atoms'], p['eqs'], p['max_killable'], p['floor_failing'],
                 p['implied_score']), flush=True)
        json.dump(out, open('runs/price.json', 'w'), indent=1)
