"""S11 step 79: kernel directions + EXACT univariate polynomials along them.

unipoly.py settled the single-coordinate question exactly, with no linearisation:
along every free input that reaches a failing check, the gcd of all the checks that
currently hold has degree 1 -- so d = 0 is the ONLY jump that breaks nothing.
Single-coordinate freedom does not exist.  The freedom, if any, is multi-coordinate.

Interpolating a general multivariate polynomial is hopeless, but we do not need one.
Use the two tools together:

  * the mod-p Jacobian picks the DIRECTION -- a kernel vector w of the rows that must
    stay zero, so the first-order term along w vanishes identically;
  * the exact univariate interpolation along that single direction then says what
    actually happens at finite distance, where every ceiling in this lab is silent.

Along a kernel direction f_c(d) has a double root at 0 for every check that holds.
gcd over those checks, divided by x^2, has as its roots the FINITE jumps that break
nothing -- nonlinear symmetries of the solution set.  Intersecting with the failing
checks' polynomials gives the jumps that fix them.

Usage: kerpoly.py START END [state.json] [K] [NDIR]
"""
import os, sys, time, random, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
import suppfree
from chunk import sweep, load
import fpoly as UP
P = ad.P
random.seed(5)

src = sys.argv[3] if len(sys.argv) > 3 else 'mod9118_0.json'
K = int(sys.argv[4]) if len(sys.argv) > 4 else 9
NDIR = int(sys.argv[5]) if len(sys.argv) > 5 else 200
tag = 'ker_' + os.path.basename(src).replace('.json', '')
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
bm = [x % P for x in base]
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
CHECKS = [a for a in range(L.NA) if a not in L.atom_out]
BAD = [a for a in CHECKS if bav[a]]
print(f'{src}: score {BASE}; failing checks {BAD}', flush=True)

_, freelist, SVS = suppfree.build(base, modp=None)
U = set()
for c in BAD[:2]:
    m = suppfree.atom_supp(c, base, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
print(f'{len(U)} columns (structural support of {BAD[:2]})', flush=True)

t0 = time.time()
cols = {u: jac_column(u, base, bm, CHECKS) for u in U}
rows = sorted(set().union(*[set(c) for c in cols.values()]))
ZERO = [c for c in rows if bav[c] == 0]
print(f'jacobian {len(rows)} rows x {len(U)} cols ({time.time()-t0:.0f}s); '
      f'{len(ZERO)} of them currently hold', flush=True)

# kernel of the must-stay-zero rows
m = len(U)
A = [[cols[u].get(c, 0) % P for u in U] for c in ZERO]
piv, r_ = [], 0
n = len(A)
for j in range(m):
    k = next((i for i in range(r_, n) if A[i][j]), None)
    if k is None:
        continue
    A[r_], A[k] = A[k], A[r_]
    inv = pow(A[r_][j], -1, P)
    A[r_] = [x * inv % P for x in A[r_]]
    for i in range(n):
        if i != r_ and A[i][j]:
            f = A[i][j]
            A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
    piv.append(j)
    r_ += 1
free_cols = [j for j in range(m) if j not in set(piv)]
KER = []
for j0 in free_cols:
    w = [0] * m
    w[j0] = 1
    for i, j in enumerate(piv):
        w[j] = (-A[i][j0]) % P
    KER.append(w)
print(f'rank {r_}; KERNEL DIMENSION {len(KER)} '
      f'(directions whose first-order effect on every holding check is zero)',
      flush=True)

DIRS = [('basis%d' % i, w) for i, w in enumerate(KER)]
for i in range(NDIR - len(KER)):
    if not KER:
        break
    w = [0] * m
    for b in random.sample(KER, min(len(KER), random.randint(2, 5))):
        s = random.randrange(1, P)
        w = [(x + s * y) % P for x, y in zip(w, b)]
    DIRS.append(('rand%d' % i, w))
print(f'{len(DIRS)} directions to scan', flush=True)


def evaluate(spec):
    name, w = spec
    xs = list(range(K + 3))
    cols_v = []
    for d in xs:
        v = list(base)
        for j, u in enumerate(U):
            if w[j]:
                v[u] = v[u] + d * w[j]
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        cols_v.append([av[c] % P for c in CHECKS])
    vary = [i for i in range(len(CHECKS)) if len(set(cv[i] for cv in cols_v)) > 1]
    G, fails, degbad = None, [], 0
    for i in vary:
        f = UP.interp(xs[:K + 1], [cols_v[d][i] for d in xs[:K + 1]])
        if not all(sum(co * pow(x, e, P) for e, co in enumerate(f)) % P == cols_v[x][i]
                   for x in xs[K + 1:]):
            degbad += 1
            continue
        c = CHECKS[i]
        if bav[c] == 0:
            G = f if G is None else UP.pgcd(G, f)
        else:
            fails.append(f)
    out = {'name': name, 'vary': len(vary), 'degfail': degbad,
           'Gdeg': (len(G) - 1) if G else None}
    H = G
    for f in fails:
        H = f if H is None else UP.pgcd(H, f)
    for label, poly in (('safe', G), ('fix', H)):
        if poly and len(poly) > 1:
            rs = [r for r in UP.roots(poly) if r % P]
            out[label + 'roots'] = len(rs)
            if rs and label == 'fix':
                best = BASE
                for r in rs[:6]:
                    v = list(base)
                    for j, u in enumerate(U):
                        if w[j]:
                            v[u] = v[u] + r * w[j]
                    ad.fwd(v, rounds=6)
                    av = L.all_atom_values(v)
                    s = L.NEQ - len(L.failing_eqs(av))
                    if s > best:
                        best = s
                        T.save(v, os.path.join(HERE, 'KER_%d_%s.json' % (s, name)))
                out['score'] = best
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(DIRS)
sweep(tag, DIRS, evaluate, start, min(end, len(DIRS)),
      keyfn=lambda s: s[0], budget=540)
rs = load(tag)
if rs:
    print('\ndirections with a finite jump that breaks nothing:')
    for r in sorted(rs, key=lambda r: -(r.get('saferoots') or 0))[:15]:
        print('   %-10s varying %-4d gcd deg %-4s safe-roots %-4s fix-roots %-4s '
              'score %s' % (r['name'], r['vary'], r.get('Gdeg'),
                            r.get('saferoots'), r.get('fixroots'), r.get('score')))
