"""S10 step 66: is 7 an INVARIANT rather than a property of the placement?

eighth.py found: adding an 8th adjustable atom enlarges the region 12 -> 15 (or 17)
equations and raises the max satisfied 5 -> 8 (or 10) -- failing stays EXACTLY 7.
The extra free parameter buys precisely as many equations as it drags in.

Test that with two and three extra atoms, and with the largest region reachable.
"""
import os, sys, collections, math, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
M1 = 7376877 * P
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E0 = sorted(L.eqs_of_atoms(NZ))
D0 = v[7068] - v[2099]
K2 = v[28730] % P

SOLO = collections.defaultdict(list)
for u in range(L.NVARS):
    if u not in L.definer and len(L.var_atoms[u]) == 1:
        SOLO[L.var_atoms[u][0]].append(u)
GRAN = {}
for a, us in SOLO.items():
    g = 0
    for u in us:
        r = T.lin_parts(a, u, v)
        if r:
            g = math.gcd(g, abs(r[0]))
    if g:
        GRAN[a] = g


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


def solve_int(rows_, rhs, n):
    m = len(rows_)
    A = [r[:] for r in rows_]
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
        piv.append(nz[0] if nz else None)
    w = [0] * n; b = list(rhs)
    for r in range(m):
        c = piv[r]
        if c is None:
            if b[r] != 0: return None
            continue
        if b[r] % A[r][c]: return None
        w[c] = b[r] // A[r][c]
        for rr in range(r + 1, m): b[rr] -= A[rr][c] * w[c]
    return [sum(U[i][c] * w[c] for c in range(n)) for i in range(n)]


def best_for(extra):
    ATOMS = NZ + list(extra)
    EQS = sorted(set(E0) | set().union(*[set(L.atom2eq.get(x, ())) for x in extra])
                 if extra else set(E0))
    M = [[L.eq_atoms[e][2].get(a, 0) for a in ATOMS] for e in EQS]
    n = len(ATOMS)
    for size in range(len(EQS), 0, -1):
        for S in itertools.combinations(range(len(EQS)), size):
            B = int_kernel([M[i] for i in S])
            if not B:
                continue
            rows_ = [[b[0] + 7376877 * b[6] for b in B] + [M1] + [0] * (1 + len(extra)),
                     [b[1] for b in B] + [0, P] + [0] * len(extra)]
            rhs = [D0, K2]
            for t, X in enumerate(extra):
                r = [b[7 + t] for b in B] + [0, 0] + [0] * len(extra)
                r[len(B) + 2 + t] = GRAN[X]
                rows_.append(r); rhs.append(0)
            if solve_int(rows_, rhs, len(B) + 2 + len(extra)) is not None:
                return len(EQS), size, len(EQS) - size
        if size < len(EQS) - 16:
            break
    return len(EQS), 0, len(EQS)


CAND = [35756, 35754]
print(f'{"extra atoms":>24} {"region":>7} {"satisfied":>10} {"FAILING":>8}')
for extra in ([], [35756], [35754], [35756, 35754]):
    R, S, F = best_for(extra)
    print(f'{str(extra):>24} {R:>7} {S:>10} {F:>8}', flush=True)
