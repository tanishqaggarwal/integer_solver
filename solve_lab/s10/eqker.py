"""S11 step 80: the kernel of the EQUATIONS, not of the atoms.

kerpoly found a 43-dimensional kernel and every direction in it moved nothing at
all: the kernel of the atom-level Jacobian is exactly the span of the handle
columns, which are identically zero mod p.  On the real inputs the atom-level
Jacobian has FULL rank -- 154 of 154 -- so at the atom level the state is rigid.

But atoms do not have to be zero.  An equation is satisfied iff its atom
COMBINATION vanishes, and every closure in this lab has been built from atom rows,
over-constraining the problem by the whole difference.  Two things change with
equation rows:

  * a satisfied equation contributes ONE row (its combination) instead of one row
    per atom -- 39,033 possible rows instead of 42,267, and far fewer per closure;
  * a SQUARED equation m*(sum)^2 has derivative 2m*(sum)*d(sum), and at a satisfied
    point sum = 0, so it contributes NO first-order constraint at all.

So the equation-level kernel is strictly larger, and the exact univariate scan along
its directions is the honest question.  Rows here are the combinations S_e = sum_a
c_a * a(x); the equation holds iff S_e = 0, for the linear and squared forms alike,
so interpolating S_e keeps the degree down as well.

Usage: eqker.py START END [state.json] [K] [NDIR]
"""
import os, sys, time, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
import suppfree
import fpoly as UP
from chunk import sweep, load
P = ad.P
random.seed(17)

src = sys.argv[3] if len(sys.argv) > 3 else 'mod9118_0.json'
K = int(sys.argv[4]) if len(sys.argv) > 4 else 9
NDIR = int(sys.argv[5]) if len(sys.argv) > 5 else 250
tag = 'eqker_' + os.path.basename(src).replace('.json', '')
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
bm = [x % P for x in base]
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
CHECKS = [a for a in range(L.NA) if a not in L.atom_out]
BAD = [a for a in CHECKS if bav[a]]
FAILEQ = sorted(L.failing_eqs(bav))
print(f'{src}: score {BASE}; failing checks {BAD}; failing equations {len(FAILEQ)}',
      flush=True)

_, freelist, SVS = suppfree.build(base, modp=None)
U = set()
for c in BAD:
    m = suppfree.atom_supp(c, base, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
t0 = time.time()
cols = {u: jac_column(u, base, bm, CHECKS) for u in U}
U = [u for u in U if cols[u]]                      # drop the handle columns
print(f'{len(U)} non-handle columns ({time.time()-t0:.0f}s)', flush=True)

atoms_moved = sorted(set().union(*[set(cols[u]) for u in U]))
EQS = sorted(set().union(*[set(L.atom2eq[a]) for a in atoms_moved]))
HOLD = [e for e in EQS if e not in set(FAILEQ)]
print(f'{len(atoms_moved)} atoms move -> {len(EQS)} equations touched, '
      f'{len(HOLD)} of them currently hold', flush=True)


def eqrow(e, u):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        d = cols[u].get(a, 0)
        if d:
            s += c * d
    return s % P


m = len(U)
A = [[eqrow(e, u) for u in U] for e in HOLD]
A = [r for r in A if any(r)]
print(f'equation-level jacobian: {len(A)} nonzero rows x {m} cols', flush=True)
piv, r_, n = [], 0, len(A)
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
KER = []
for j0 in [j for j in range(m) if j not in set(piv)]:
    w = [0] * m
    w[j0] = 1
    for i, j in enumerate(piv):
        w[j] = (-A[i][j0]) % P
    KER.append(w)
print(f'rank {r_};  EQUATION-LEVEL KERNEL DIMENSION {len(KER)}  '
      f'(atom-level kernel on the same columns was 0)', flush=True)

DIRS = [('basis%d' % i, w) for i, w in enumerate(KER)]
for i in range(max(0, NDIR - len(KER))):
    if not KER:
        break
    w = [0] * m
    for b in random.sample(KER, min(len(KER), random.randint(2, 6))):
        s = random.randrange(1, P)
        w = [(x + s * y) % P for x, y in zip(w, b)]
    DIRS.append(('rand%d' % i, w))
print(f'{len(DIRS)} directions to scan exactly', flush=True)
ALLEQ = sorted(set(EQS) | set(FAILEQ))


def comb(av, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if av[a]:
            s += c * av[a]
    return s % P


def evaluate(spec):
    name, w = spec
    xs = list(range(K + 3))
    vals = []
    for d in xs:
        v = list(base)
        for j, u in enumerate(U):
            if w[j]:
                v[u] = v[u] + d * w[j]
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        vals.append([comb(av, e) for e in ALLEQ])
    vary = [i for i in range(len(ALLEQ)) if len(set(x[i] for x in vals)) > 1]
    G, fails, degbad = None, [], 0
    for i in vary:
        f = UP.interp(xs[:K + 1], [vals[d][i] for d in xs[:K + 1]])
        if not all(sum(co * pow(x, e2, P) for e2, co in enumerate(f)) % P == vals[x][i]
                   for x in xs[K + 1:]):
            degbad += 1
            continue
        if vals[0][i] == 0:
            G = f if G is None else UP.pgcd(G, f)
        else:
            fails.append(f)
    out = {'name': name, 'vary': len(vary), 'degfail': degbad,
           'Gdeg': (len(G) - 1) if G else None, 'nfail': len(fails)}
    H = G
    for f in fails:
        H = f if H is None else UP.pgcd(H, f)
    for label, poly in (('safe', G), ('fix', H)):
        if poly and len(poly) > 1:
            rs = [r for r in UP.roots(poly) if r % P]
            out[label + 'roots'] = len(rs)
            if rs:
                best = BASE
                for r in rs[:6]:
                    v = list(base)
                    for j, u in enumerate(U):
                        if w[j]:
                            v[u] = v[u] + r * w[j]
                    ad.fwd(v, rounds=6)
                    s = L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))
                    if s > best:
                        best = s
                        T.save(v, os.path.join(HERE, 'EQK_%d_%s.json' % (s, name)))
                out[label + 'score'] = best
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(DIRS)
sweep(tag, DIRS, evaluate, start, min(end, len(DIRS)), keyfn=lambda s: s[0],
      budget=540)
rs = load(tag)
if rs:
    print('\nbest directions:')
    for r in sorted(rs, key=lambda r: -max(r.get('safescore') or 0,
                                           r.get('fixscore') or 0))[:15]:
        print('   %-10s vary %-4d Gdeg %-4s safe-roots %-4s fix-roots %-4s '
              'score %s/%s' % (r['name'], r['vary'], r.get('Gdeg'),
                               r.get('saferoots'), r.get('fixroots'),
                               r.get('safescore'), r.get('fixscore')))
