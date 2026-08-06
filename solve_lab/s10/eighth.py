"""S10 step 64: enlarge the defect placement by one atom.

Part I proved 5-of-12 optimal for the SEVEN-atom placement with its two mod-p
congruences.  That proof is conditional on the placement.  Adding an 8th
adjustable atom X adds one free parameter -- worth up to +1 satisfied equation --
at the price of the equations X drags in from outside the 12.

So: find adjustable atoms whose equation footprint is (almost) inside the 12.
Any with footprint fully inside is a FREE extra parameter -> 39,027 immediately.
"""
import os, sys, collections, math, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = L.load(BEST)
av = L.all_atom_values(v)
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = sorted(L.eqs_of_atoms(NZ))
ES = set(E)
print(f'current placement: {len(NZ)} atoms, {len(E)} equations, 7 failing')

# adjustable atoms: carry a free solo handle
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
print(f'adjustable atoms: {len(GRAN)}')

# ALSO: any atom that is currently zero and whose variables include a free input
# occurring nowhere else that matters
rows = []
for a, g in GRAN.items():
    if a in NZ or av[a]:
        continue
    eqs = set(L.atom2eq.get(a, ()))
    if not eqs:
        continue
    over = eqs - ES
    inside = eqs & ES
    if inside:
        rows.append((len(over), len(inside), a, g))
rows.sort()
print(f'\nadjustable atoms touching the 12 equations: {len(rows)}')
print(f'{"overflow":>9} {"inside":>7} {"atom":>7} {"granularity":>13}')
for over, ins, a, g in rows[:25]:
    gg = f'{g//P}*p' if g % P == 0 else str(g)
    print(f'{over:>9} {ins:>7} {a:>7} {gg:>13}   {L.atom_src[a][:60]}')

if rows and rows[0][0] == 0:
    print('\n*** an adjustable atom lies ENTIRELY inside the 12 equations -> free parameter')
else:
    print(f'\nminimum overflow = {rows[0][0] if rows else "-"} extra equations')

# ---- exact re-optimisation with the best few candidates added ---------------
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


D0 = v[7068] - v[2099]
K2 = v[28730] % P
M1 = 7376877 * P
print('\n=== exact re-optimisation with one extra atom ===')
for over, ins, X, g in rows[:8]:
    ATOMS = NZ + [X]
    EQS = sorted(ES | set(L.atom2eq.get(X, ())))
    M = [[L.eq_atoms[e][2].get(a, 0) for a in ATOMS] for e in EQS]
    n = len(ATOMS)
    bestsz = 0
    for size in range(len(EQS), 0, -1):
        ok = False
        for S in itertools.combinations(range(len(EQS)), size):
            sub = [M[i] for i in S]
            B = int_kernel(sub)
            if not B: continue
            c1 = [b[0] + 7376877 * b[6] for b in B]
            c2 = [b[1] for b in B]
            c3 = [b[7] for b in B]          # extra atom must be a multiple of g
            sol = solve_int([c1 + [M1, 0, 0], c2 + [0, P, 0], c3 + [0, 0, g]],
                            [D0, K2, 0], len(B) + 3)
            if sol is not None:
                ok = True; break
        if ok:
            bestsz = size; break
        if size < len(EQS) - 14:
            break
    print(f'  atom {X:<7} overflow={over} region={len(EQS)} eqs -> max satisfied '
          f'{bestsz}  => failing {len(EQS)-bestsz}  (current 7)', flush=True)
