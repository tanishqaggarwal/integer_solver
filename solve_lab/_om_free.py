#!/usr/bin/env python3
"""Which atoms are FORCED to zero by the equations (over Q, tested mod a prime)?
Builds the closure around the failing region, dense-eliminates mod P, and reports
the atoms that can be nonzero in some solution supported on the closure."""
import pickle, sys
import numpy as np
from collections import defaultdict

P = 2147483647  # 2^31-1

D = pickle.load(open('_om_parsed2.pkl', 'rb'))
eqatoms = D['eqatoms']
ainc = defaultdict(list)
for e, d in enumerate(eqatoms):
    for k in d: ainc[k].append(e)
F = [2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
ALPHA = '((x7068-x2099)-(7376877*x642))'
BETA = '((x4432-x19964)-x28730)'

def closure(shells):
    A = set()
    for e in F: A |= set(eqatoms[e])
    for s in range(shells):
        E = set(F)
        for k in A: E |= set(ainc[k])
        newA = set(A)
        for e in E: newA |= set(eqatoms[e])
        if newA == A:
            print('closed at shell', s); break
        A = newA
    E = set(F)
    for k in A: E |= set(ainc[k])
    return sorted(A), sorted(E)

def rref_mod(M, P):
    """in-place RREF mod P; returns pivot column list"""
    nr, nc = M.shape
    piv = []; r = 0
    for c in range(nc):
        if r >= nr: break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), P - 2, P)
        M[r] = (M[r] * inv) % P
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % P
        piv.append(c); r += 1
    return piv

if __name__ == '__main__':
    shells = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    A, E = closure(shells)
    print('|A|=%d |E|=%d' % (len(A), len(E)))
    idx = {k: i for i, k in enumerate(A)}
    M = np.zeros((len(E), len(A)), dtype=np.int64)
    for r, e in enumerate(E):
        for k, c in eqatoms[e].items():
            if k in idx: M[r, idx[k]] = c % P
    piv = rref_mod(M, P)
    rk = len(piv)
    print('rank=%d nullity=%d' % (rk, len(A) - rk))
    pivset = set(piv)
    free = [c for c in range(len(A)) if c not in pivset]
    # nullspace: for each free col f, vector with 1 at f and -M[i,f] at piv[i]
    # atom a can be nonzero  <=>  a is free, OR some free col f has M[row_of_a, f] != 0
    canbe = set(free)
    for i, c in enumerate(piv):
        row = M[i, free]
        if np.any(row != 0): canbe.add(c)
    print('atoms that CAN be nonzero: %d / %d' % (len(canbe), len(A)))
    print('alpha can be nonzero:', idx[ALPHA] in canbe)
    print('beta  can be nonzero:', idx[BETA] in canbe)
    print()
    print('sample of atoms that can be nonzero:')
    for c in sorted(canbe)[:40]:
        print('   ', A[c])
