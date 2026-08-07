"""S11 step 114: Part XXXI's machinery, applied to the ALGEBRAIC states.

Part XXXI settled the coding optimum exactly -- at the deliverable the twelve
equations are homogeneous in the seven residual atom values, the zero-cost generators
are enumerable (nine of fifteen candidates), and the Diophantine problem over that
lattice caps the satisfiable rows at five, giving 39,026.

The same three steps apply to any state.  Run them on the algebraic states, where the
the two-condition primitive closes and the residual is a different, smaller object:

  1. residual atoms R, and the equations E they touch;  is E homogeneous in R?
  2. every free input reaching R: its exact effect on each atom of R, and its
     MEASURED cost (equations that break outside E) -- never estimated from atom
     counts, which was wrong twice;
  3. for every subset of E, the linear Diophantine system over the zero-cost lattice,
     largest subsets first.

Usage: algside.py [state.json] [MAXCAND]
"""
import os, sys, time, itertools, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'AG_39013.json'
MAXC = int(sys.argv[2]) if len(sys.argv) > 2 else 120
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
av0 = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(av0))
R = [a for a in range(L.NA) if a not in L.atom_out and av0[a]]
E = sorted(set(e for a in R for e in L.atom2eq[a]))
FAILE = sorted(L.failing_eqs(av0))
print('%s: score %d; residual atoms %s' % (src, BASE, R), flush=True)
print('   they touch %d equations, of which %d fail'
      % (len(E), len([e for e in E if e in set(FAILE)])), flush=True)
print('   failing equations not touched by them: %d'
      % len(set(FAILE) - set(E)), flush=True)

# homogeneity: are the other atoms of those equations zero?
extra = sorted({a for e in E for a in L.eq_atoms[e][2] if a not in set(R)
                and av0[a]})
print('   nonzero atoms in E outside R: %s' % extra, flush=True)
RR = sorted(set(R) | set(extra))
M = [[L.eq_atoms[e][2].get(a, 0) for a in RR] for e in E]
alpha0 = [av0[a] for a in RR]
NSAT = sum(1 for i in range(len(E))
           if sum(M[i][j] * alpha0[j] for j in range(len(RR))) == 0)
print('   system: %d rows x %d atom values; currently satisfied %d'
      % (len(E), len(RR), NSAT), flush=True)

# candidate generators: free inputs reaching any atom of RR
_, freelist, SVS = suppfree.build(base, modp=None)
U = set()
for a in RR:
    m = suppfree.atom_supp(a, base, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
print('   %d free inputs reach them; measuring the cheapest %d'
      % (len(U), min(MAXC, len(U))), flush=True)
ES = set(E)
GEN, t0 = [], time.time()
for i, u in enumerate(U[:MAXC]):
    v = list(base)
    v[u] = v[u] + 1
    ad.fwd(v, rounds=6)
    aw = L.all_atom_values(v)
    d = [aw[a] - av0[a] for a in RR]
    cost = len(set(L.failing_eqs(aw)) - ES)
    if cost == 0 and any(d):
        GEN.append((u, d))
    if i % 25 == 0:
        print('      %d/%d  (%.0fs)  free so far %d'
              % (i, min(MAXC, len(U)), time.time() - t0, len(GEN)), flush=True)
print('\n%d ZERO-COST generators:' % len(GEN), flush=True)
NM = ['a%d' % a for a in RR]
for u, d in GEN:
    print('   x%-6d %s' % (u, ', '.join(
        '%s%+s' % (NM[i], str(x)[:16]) for i, x in enumerate(d) if x)), flush=True)
json.dump({'state': src, 'R': R, 'E': E, 'RR': RR,
           'gen': [[u, [str(x) for x in d]] for u, d in GEN]},
          open(os.path.join(HERE, 'algside_%s.json'
                            % os.path.basename(src).replace('.json', '')), 'w'))
if not GEN:
    print('no zero-cost generator: nothing to solve'); sys.exit()

A = [[sum(M[i][j] * d[j] for j in range(len(RR))) for _, d in GEN]
     for i in range(len(E))]
B = [-sum(M[i][j] * alpha0[j] for j in range(len(RR))) for i in range(len(E))]


def int_solve(rows, rhs):
    n, m = len(rows), len(rows[0])
    H = [r[:] for r in rows]
    Umat = [[1 if i == j else 0 for j in range(m)] for i in range(m)]
    piv, r = [], 0
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
                        Umat[k][j] -= q * Umat[k][j0]
        nz = [j for j in range(r, m) if H[i][j]]
        if not nz:
            continue
        j0 = nz[0]
        if j0 != r:
            for k in range(n):
                H[k][r], H[k][j0] = H[k][j0], H[k][r]
            for k in range(m):
                Umat[k][r], Umat[k][j0] = Umat[k][j0], Umat[k][r]
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
    return [sum(Umat[k][j] * y[j] for j in range(m)) for k in range(m)]


print('\nDiophantine search over %d rows and %d generators' % (len(E), len(GEN)),
      flush=True)
t0 = time.time()
for k in range(len(E), NSAT, -1):
    if time.time() - t0 > 900:
        print('   (time budget reached at |S| = %d)' % k, flush=True)
        break
    hits = 0
    for S in itertools.combinations(range(len(E)), k):
        c = int_solve([A[i] for i in S], [B[i] for i in S])
        if c is None:
            continue
        hits += 1
        v = list(base)
        for (u, _), ci in zip(GEN, c):
            if ci:
                v[u] = v[u] + ci
        ad.fwd(v, rounds=6)
        aw = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(aw))
        print('   |S|=%d -> SCORE %d  rows %s'
              % (k, s, [E[i] for i in S][:8]), flush=True)
        if s > BASE:
            T.save(v, os.path.join(HERE, 'ALG_%d.json' % s))
            print('      *** NEW BEST %d -- saved ALG_%d.json' % (s, s), flush=True)
        if hits >= 8:
            break
    print('|S| = %d : %d integrally solvable  (%.0fs)'
          % (k, hits, time.time() - t0), flush=True)
    if hits:
        break
