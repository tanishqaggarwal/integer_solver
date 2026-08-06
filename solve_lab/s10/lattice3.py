"""S10 step 18: residual optimum with the p-shift of x_7068 included.

repairD.py verified: x_7068 += k*p, then repair a29539 via x_29967 and a40826 via
x_30163, leaves EXACTLY the same 7 nonzero atoms and nothing else.  So D is free
within D0 + p*Z, and the constraint

    A1 + 7376877*A7 == D  (mod 7376877*p)   for SOME D in D0 + p*Z

collapses to the single weaker congruence

    A1 + 7376877*A7 == D0   (mod p)          <-- 7376877 factor gone

together with   A2 == K2 (mod p).
"""
import os, sys, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = L.load(BEST)
av = L.all_atom_values(v)
EQS = sorted(L.eqs_of_atoms(NZ))
M = [[L.eq_atoms[i][2].get(a, 0) for a in NZ] for i in EQS]
D0 = v[7068] - v[2099]
K2 = v[28730] % P
print('D0 % p =', D0 % P)
print('K2     =', K2)


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
    w = [0] * nvars; b = list(rhs)
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
    sol = solve_int_system([c1 + [P, 0], c2 + [0, P]], [D0 % P, K2], k + 2)
    if sol is None: return None
    y = sol[:k]
    A = [sum(y[j] * B[j][i] for j in range(k)) for i in range(7)]
    for i in S:
        if sum(M[EQS.index(i)][t] * A[t] for t in range(7)) != 0: return None
    if (A[0] + 7376877 * A[6] - D0) % P or (A[1] - K2) % P: return None
    return A


best = None
for size in range(12, 0, -1):
    found = []
    for S in itertools.combinations(EQS, size):
        A = witness(list(S))
        if A is not None:
            found.append((S, A))
            if len(found) >= 300: break
    print(f'  size {size:>2}: {len(found)} integer-solvable')
    if found:
        best = (size, found); break

size, found = best
print(f'\nOPTIMUM with the p-shift: {size} of {len(EQS)} -> score {L.NEQ-len(EQS)+size}')
json.dump({'size': size, 'EQS': EQS, 'NZ': NZ, 'D0': str(D0), 'K2': str(K2),
           'cands': [{'S': list(S), 'A': [str(x) for x in A]} for S, A in found]},
          open(os.path.join(HERE, 'lattice3_best.json'), 'w'))
print('saved lattice3_best.json')
