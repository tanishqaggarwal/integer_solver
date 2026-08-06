"""S10 step 68: a22231 is a FREE 8th atom -- a gap in my own model.

isolate.py moved x_28730 together with x_4432 so that a22231 stayed 0, and THAT
is what dragged in atom 7930 (15 equations).  But a22231 need not be zero.  Move
x_28730 ALONE:

    a22230 = x_28730 - p*x_9413                 changes by +d
    a22231 = x_4432 - x_19964 - x_28730         changes by -d
    x_4432 untouched  ->  NO downstream collateral

So with K = x_4432 - x_19964 fixed:
    A2 = a22230 = t8 - p*t3       A8 = a22231 = K - t8
    => A2 + A8 == K (mod p)   with A2 - A8 otherwise FREE

That is ONE congruence on EIGHT atoms instead of one on seven -- strictly more
freedom than the model Part I optimised.
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
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = L.load(BEST)
av = L.all_atom_values(v)
base_nz = set(a for a in range(L.NA) if av[a])

print('a22231 =', L.atom_src[22231])
print('  currently', av[22231], ' in', len(L.atom2eq.get(22231, {})), 'equations:',
      sorted(L.atom2eq.get(22231, {})))
print('a22230 =', L.atom_src[22230])
print('x_28730 atoms:', sorted(L.var_atoms[28730]))

# ---- empirical: move x_28730 ALONE, blocking a22231 from repair -------------
print('\n=== empirical: move x_28730 alone (a22231 NOT repaired) ===')
NZ7 = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ7) | {22231}
for d in (1, 2, P):
    w = list(v)
    L.ripple(w, {28730: v[28730] + d}, block=BLOCK)
    aw = L.all_atom_values(w)
    nz = set(a for a in range(L.NA) if aw[a])
    extra = sorted(nz - base_nz)
    print(f'  d={"p" if d == P else d}: new nonzero atoms {extra}  '
          f'x_4432 moved: {w[4432] != v[4432]}  failing={len(L.failing_eqs(aw))}')

# ---- exact lattice optimisation with the 8-atom model ----------------------
NZ = NZ7 + [22231]
E = sorted(L.eqs_of_atoms(NZ))
print(f'\nregion: {len(NZ)} atoms, {len(E)} equations  {E}')
M = [[L.eq_atoms[e][2].get(a, 0) for a in NZ] for e in E]
D0 = v[7068] - v[2099]
K = v[4432] - v[19964]
print(f'K = x_4432 - x_19964 ; K mod p = {K % P}')
print(f'sanity A2+A8-K mod p = {(av[22230] + av[22231] - K) % P}')


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


# constraints:  A1 + 7376877*A7 == D0 (mod p)   [p-shift of x_7068 is free]
#               A2 + A8         == K  (mod p)
print('\n=== exact optimisation, 8-atom model ===')
found = None
for size in range(len(E), 0, -1):
    hit = None
    for S in itertools.combinations(range(len(E)), size):
        B = int_kernel([M[i] for i in S])
        if not B:
            continue
        c1 = [b[0] + 7376877 * b[6] for b in B]
        c2 = [b[1] + b[7] for b in B]
        sol = solve_int([c1 + [P, 0], c2 + [0, P]], [D0 % P, K % P], len(B) + 2)
        if sol is not None:
            y = sol[:len(B)]
            A = [sum(y[j] * B[j][i] for j in range(len(B))) for i in range(len(NZ))]
            if any(A):
                hit = (S, A); break
    print(f'  size {size:>2}: {"SOLVABLE" if hit else "no"}', flush=True)
    if hit:
        found = (size, hit); break

if found:
    size, (S, A) = found
    print(f'\n*** max satisfied {size} of {len(E)} -> failing {len(E)-size} '
          f'(score {L.NEQ-(len(E)-size)})')
    print(f'    satisfied equations: {[E[i] for i in S]}')
    print(f'    atom values: {[str(x)[:26] for x in A]}')
    json.dump({'NZ': NZ, 'E': E, 'S': [E[i] for i in S], 'A': [str(x) for x in A],
               'K': str(K), 'D0': str(D0)},
              open(os.path.join(HERE, 'a22231_best.json'), 'w'))
    print('    saved a22231_best.json')
