"""S10 step 15: EXACT residual optimum with the CORRECT achievable set.

Verified knob freedom (isolate.py): x_642, x_17325, x_9413, x_1329, x_10903,
x_29854, x_31864, x_9118, x_8731 are free (no collateral atoms).  x_28730 is NOT
(it drags x_4432 and breaks atoms 7930 + 41512).

=> achievable set of A = (a22229,a22230,a35758,a35759,a35760,a35761,a35762):
     (1)  A1 + 7376877*A7 == D    (mod 7376877*p)
     (2)  A2             == K2    (mod p)          [x_28730 pinned, x_9413 free]
     (3)  A3,A4,A5,A6 free
"""
import os, sys, json, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
M1 = 7376877 * P
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = L.load(BEST)
av = L.all_atom_values(v)
EQS = sorted(L.eqs_of_atoms(NZ))
M = [[L.eq_atoms[i][2].get(a, 0) for a in NZ] for i in EQS]
D = v[7068] - v[2099]
K2 = v[28730] % P
A0 = [av[a] for a in NZ]
print('sanity 1:', (A0[0] + 7376877 * A0[6] - D) % M1)
print('sanity 2:', (A0[1] - K2) % P)
sat0 = [i for i in EQS if L.eq_value(i, av) == 0]
print(f'base satisfies {len(sat0)}/{len(EQS)} -> score {L.NEQ-len(EQS)+len(sat0)}')


def int_kernel(mat):
    m, n = len(mat), len(mat[0])
    A = [r[:] for r in mat]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    piv = []
    for r in range(m):
        while True:
            nz = [c for c in range(n) if c not in piv and A[r][c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda c: abs(A[r][c])); p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m): A[i][c] -= q * A[i][p0]
                    for i in range(n): U[i][c] -= q * U[i][p0]
        nz = [c for c in range(n) if c not in piv and A[r][c] != 0]
        if nz: piv.append(nz[0])
    return [[U[i][c] for i in range(n)] for c in range(n) if c not in piv]


def solve_int_system(rows, rhs, nvars):
    """Integer solution x of rows.x = rhs (2 rows), or None.  Column HNF."""
    m = len(rows)
    A = [r[:] for r in rows]
    U = [[1 if i == j else 0 for j in range(nvars)] for i in range(nvars)]
    piv = []
    for r in range(m):
        while True:
            nz = [c for c in range(nvars) if c not in piv and A[r][c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda c: abs(A[r][c])); p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m): A[i][c] -= q * A[i][p0]
                    for i in range(nvars): U[i][c] -= q * U[i][p0]
        nz = [c for c in range(nvars) if c not in piv and A[r][c] != 0]
        piv.append(nz[0] if nz else None)
    w = [0] * nvars
    b = list(rhs)
    for r in range(m):
        c = piv[r]
        if c is None:
            if b[r] != 0: return None
            continue
        if b[r] % A[r][c]: return None
        w[c] = b[r] // A[r][c]
        for rr in range(r + 1, m): b[rr] -= A[rr][c] * w[c]
    return [sum(U[i][c] * w[c] for c in range(nvars)) for i in range(nvars)]


def witness(S):
    rows = [M[EQS.index(i)] for i in S]
    B = int_kernel(rows) if rows else [[1 if i == j else 0 for i in range(7)] for j in range(7)]
    if not B: return None
    k = len(B)
    c1 = [b[0] + 7376877 * b[6] for b in B]
    c2 = [b[1] for b in B]
    sol = solve_int_system([c1 + [M1, 0], c2 + [0, P]], [D, K2], k + 2)
    if sol is None: return None
    y = sol[:k]
    A = [sum(y[j] * B[j][i] for j in range(k)) for i in range(7)]
    for i in S:
        if sum(M[EQS.index(i)][t] * A[t] for t in range(7)) != 0: return None
    if (A[0] + 7376877 * A[6] - D) % M1 or (A[1] - K2) % P: return None
    return A


best = None
for size in range(12, 0, -1):
    found = []
    for S in itertools.combinations(EQS, size):
        A = witness(list(S))
        if A is not None:
            found.append((S, A))
            if len(found) > 400: break
    print(f'  size {size:>2}: {len(found)} integer-solvable')
    if found:
        best = (size, found)
        break

size, found = best
print(f'\nEXACT OPTIMUM for this defect placement: {size} of {len(EQS)} '
      f'-> score {L.NEQ-len(EQS)+size}')
cands = [{'S': list(S), 'A': [str(x) for x in A]} for S, A in found[:400]]
json.dump({'size': size, 'cands': cands, 'EQS': EQS, 'NZ': NZ},
          open(os.path.join(HERE, 'lattice2_best.json'), 'w'))
print(f'saved {len(cands)} candidate targets')
