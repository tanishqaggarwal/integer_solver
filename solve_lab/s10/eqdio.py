"""S11 step 107: two liftable equations, two handles -- a linear Diophantine solve.

eqlift finds the gap but cannot close it one move at a time: at PF_best_39015 two
equations (7469 and 21382) have combinations that vanish mod p but not over Z, and
fixing either with the handle x30317 breaks the other through the same handle.  With
strict improvement required, no single move exists.

The fix is to solve them TOGETHER.  For a set of knobs the requirement is

    sum_u  g[e][u] * d_u  =  -S_e     for e in {7469, 21382},   d_u integer

with g[e][u] the exact integer coefficient of knob u on equation e's combination
(intad.jacZ summed against the equation's atom coefficients).  That is a linear
Diophantine system in two rows; over pairs of knobs it is a 2x2 determinant test, so
enumerate the pairs, solve exactly, apply, and score.

Usage: eqdio.py [state.json]
"""
import os, sys, time, itertools
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PF_best_39015.json'
v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v0, rounds=6)


def comb(av, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if av[a]:
            s += c * av[a]
    return s


av = L.all_atom_values(v0)
BASE = L.NEQ - len(L.failing_eqs(av))
FAIL = sorted(L.failing_eqs(av))
LIFT = [e for e in FAIL if comb(av, e) % P == 0]
print('%s: score %d; %d failing, %d with a combination ≡ 0 mod p: %s'
      % (src, BASE, len(FAIL), len(LIFT), LIFT), flush=True)
if len(LIFT) < 1:
    sys.exit()
_, freelist, SVS = suppfree.build(v0, modp=None)
KN = set()
for e in LIFT:
    for a in L.eq_atoms[e][2]:
        m = suppfree.atom_supp(a, v0, SVS, modp=None)
        KN |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
KN = sorted(KN)
print('%d knobs reach those equations' % len(KN), flush=True)

t0 = time.time()
G = {}
for u in KN:
    atoms = sorted({a for e in LIFT for a in L.eq_atoms[e][2]})
    col = jacZ(u, v0, atoms)
    row = []
    for e in LIFT:
        s = 0
        for a, c in L.eq_atoms[e][2].items():
            if a in col:
                s += c * col[a]
        row.append(s)
    if any(row):
        G[u] = row
print('%d knobs with a nonzero exact integer effect (%.0fs)'
      % (len(G), time.time() - t0), flush=True)
S = [comb(av, e) for e in LIFT]
best, bestv = BASE, list(v0)
KS = sorted(G)
tried = 0
for u1, u2 in itertools.combinations(KS, 2):
    if time.time() - t0 > 900:
        break
    a1, b1 = G[u1][0], G[u1][1] if len(LIFT) > 1 else 0
    a2, b2 = G[u2][0], G[u2][1] if len(LIFT) > 1 else 0
    det = a1 * b2 - a2 * b1
    if det == 0:
        continue
    n1 = (-S[0]) * b2 - (-S[1]) * a2
    n2 = a1 * (-S[1]) - b1 * (-S[0])
    if n1 % det or n2 % det:
        continue
    d1, d2 = n1 // det, n2 // det
    w = list(v0)
    w[u1] += d1
    w[u2] += d2
    ad.fwd(w, rounds=6)
    aw = L.all_atom_values(w)
    s = L.NEQ - len(L.failing_eqs(aw))
    tried += 1
    if s > best:
        best, bestv = s, list(w)
        print('   x%d += %s, x%d += %s  ->  score %d'
              % (u1, str(d1)[:16], u2, str(d2)[:16], s), flush=True)
        T.save(w, os.path.join(HERE, 'ED_%d.json' % s))
print('\n%d integral pairs tried; best %d (was %d)' % (tried, best, BASE))
if best > BASE:
    T.save(bestv, os.path.join(HERE, 'ED_best_%d.json' % best))
    print('saved ED_best_%d.json' % best)
