"""S11 step 110: the zero-cost lattice of the seven residual atom values.

The twelve equations at the deliverable are HOMOGENEOUS linear forms in the seven
residual atom values alpha, with rank 7, so alpha = 0 satisfies all twelve while the
current alpha satisfies exactly five.  Which alpha are reachable at NO cost elsewhere
is decided by counting the atoms each realisation variable touches:

    x1329  x10903  x17325  x9413  x29854  x31864  x642   -- each in 1-2 atoms, ALL
                                                            among the seven: FREE
    x9118  33 extra equations   x8731  26   x7068  35   x28730  1   -- must be held

So the seven free variables generate a LATTICE of zero-cost moves in alpha-space:

    x642   += 1  ->  a6 += 1  and  a0 -= 7376877      (a22229 carries -7376877*x642)
    x29854 += 1  ->  a2 += 1,  a3 -= 1
    x31864 += 1  ->  a4 += 1,  a5 += 1
    x1329  += 1  ->  a2 -= p       x10903 += 1  ->  a4 -= p
    x9413  += 1  ->  a1 -= p       x17325 += 1  ->  a6 -= p

Every generator leaves a2+a3 and a5-a4 unchanged mod p -- exactly the condition for
the residues of x9118 and x8731, the only things with side effects, to stay put.  So
the reachable set is the coset alpha_cur + L and the question is LINEAR DIOPHANTINE:

    for a subset S of the twelve rows find integers c with  A c = B,
    A[i][j] = <M_i, g_j>,  B[i] = -<M_i, alpha_cur>

over Z, not merely over Q -- the equations must hold exactly.  |S| = 6 gives 39,027,
|S| = 7 gives 39,028, and so on up.

Usage: lattice7.py [MAXK]
"""
import os, sys, itertools, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import fwd
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
av = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(av))
E = sorted(set(e for a in SEVEN for e in L.atom2eq[a]))
M = [[L.eq_atoms[e][2].get(a, 0) for a in SEVEN] for e in E]
alpha0 = [av[a] for a in SEVEN]
NSAT = sum(1 for i in range(len(E))
           if sum(M[i][j] * alpha0[j] for j in range(7)) == 0)
print('witness %d; %d homogeneous rows; currently satisfied %d'
      % (BASE, len(E), NSAT), flush=True)

# MEASURED, not estimated: move each variable by 1 and count failures outside the
# twelve.  The atom-count estimate was wrong about x9118 (33 atoms, cost 0) and
# x8731 (26 atoms, cost 0); only x7068 (13) and x28730 (16) actually cost anything.
GEN = [(642,   [-7376877, 0, 0, 0, 0, 0, 1]),
       (29854, [0, 0, 1, -1, 0, 0, 0]),
       (31864, [0, 0, 0, 0, 1, 1, 0]),
       (1329,  [0, 0, -P, 0, 0, 0, 0]),
       (10903, [0, 0, 0, 0, -P, 0, 0]),
       (9413,  [0, -P, 0, 0, 0, 0, 0]),
       (17325, [0, 0, 0, 0, 0, 0, -P]),
       (9118,  [0, 0, 0, 5113045, 0, 0, 0]),
       (8731,  [0, 0, 0, 0, 0, 1, 0])]
G = [g for _, g in GEN]


def dot(u, w):
    return sum(a * b for a, b in zip(u, w))


A = [[dot(M[i], g) for g in G] for i in range(len(E))]
B = [-dot(M[i], alpha0) for i in range(len(E))]


def int_solve(rows, rhs):
    """Integer solution c of rows*c = rhs, via column-style Hermite reduction."""
    n, m = len(rows), len(rows[0])
    H = [r[:] for r in rows]
    U = [[1 if i == j else 0 for j in range(m)] for i in range(m)]
    piv = []                      # (row, column) pivots, in order
    r = 0
    for i in range(n):
        if r >= m:
            break
        while True:
            nz = [j for j in range(r, m) if H[i][j]]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda j: abs(H[i][j]))
            j0 = nz[0]
            for j in nz[1:]:
                q = H[i][j] // H[i][j0]
                if q:
                    for k in range(n):
                        H[k][j] -= q * H[k][j0]
                    for k in range(m):
                        U[k][j] -= q * U[k][j0]
        nz = [j for j in range(r, m) if H[i][j]]
        if not nz:
            continue
        j0 = nz[0]
        if j0 != r:
            for k in range(n):
                H[k][r], H[k][j0] = H[k][j0], H[k][r]
            for k in range(m):
                U[k][r], U[k][j0] = U[k][j0], U[k][r]
        piv.append((i, r))
        r += 1
    y = [0] * m
    for i, j in piv:
        s = rhs[i] - sum(H[i][k] * y[k] for k in range(j))
        if s % H[i][j]:
            return None
        y[j] = s // H[i][j]
    for i in range(n):
        if sum(H[i][k] * y[k] for k in range(m)) != rhs[i]:
            return None
    return [sum(U[k][j] * y[j] for j in range(m)) for k in range(m)]


def alpha_of(c):
    return [alpha0[j] + sum(c[k] * G[k][j] for k in range(7)) for j in range(7)]


def realise(c):
    v = list(base)
    for (var, _), ci in zip(GEN, c):
        if ci:
            v[var] = v[var] + ci
    fwd(v)
    return v


# sanity: c = 0 must reproduce the witness
v0 = realise([0] * 7)
print('sanity: c = 0 reproduces score %d'
      % (L.NEQ - len(L.failing_eqs(L.all_atom_values(v0)))), flush=True)

t0 = time.time()
MAXK = int(sys.argv[1]) if len(sys.argv) > 1 else 12
found_any = False
for k in range(min(MAXK, 12), NSAT, -1):
    hits = 0
    for S in itertools.combinations(range(12), k):
        c = int_solve([A[i] for i in S], [B[i] for i in S])
        if c is None:
            continue
        al = alpha_of(c)
        if not all(dot(M[i], al) == 0 for i in S):
            continue
        hits += 1
        v = realise(c)
        aw = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(aw))
        got = sum(1 for i in range(12)
                  if sum(M[i][j] * aw[SEVEN[j]] for j in range(7)) == 0)
        print('   |S|=%d rows %s -> %d of 12 satisfied, SCORE %d'
              % (k, [E[i] for i in S], got, s), flush=True)
        if s > BASE:
            T.save(v, os.path.join(HERE, 'LAT_%d.json' % s))
            print('      *** NEW BEST %d  -- saved LAT_%d.json' % (s, s),
                  flush=True)
        found_any = True
        if hits >= 12:
            break
    print('|S| = %d : %d integrally solvable subsets  (%.0fs)'
          % (k, hits, time.time() - t0), flush=True)
    if hits:
        break
if not found_any:
    print('no subset larger than %d is integrally reachable in the free lattice'
          % NSAT)
