"""S11 step 81: the equation-level kernel and the exact scan, AT THE DELIVERABLE.

eqker settled the canonical frame at mod9118_0: 3,324 holding equations against 415
non-handle columns, rank 415, EQUATION-LEVEL KERNEL DIMENSION 0.  Even after
dropping to equation rows (one row per equation instead of one per atom, and no row
at all from a satisfied squared equation) the state is first-order rigid.

But mod9118_0 is 39,009 and the deliverable is 39,026, where only seven equations
fail.  The witness lives on a different frame -- five gate atoms are deliberately
broken, so plain forward evaluation "repairs" them and destroys the state.  Redo the
whole analysis inside frame 2, where the witness is on-manifold, and ask the same
question there: is there a direction, and does any finite jump along it land back on
the satisfied set?

Usage: eqker2.py START END [K] [NDIR]
"""
import os, sys, time, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import suppfree
import fpoly as UP
from chunk import sweep, load
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score, grad, jac_column
P = ad.P
random.seed(23)
K = int(sys.argv[3]) if len(sys.argv) > 3 else 9
NDIR = int(sys.argv[4]) if len(sys.argv) > 4 else 250
tag = 'eqker2'

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
bm = [x % P for x in base]
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
BAD = [a for a in CHECKS if bav[a]]
FAILEQ = sorted(L.failing_eqs(bav))
print(f'frame 2 at the deliverable: score {BASE}; failing checks {BAD}; '
      f'failing equations {FAILEQ}', flush=True)

_, freelist, SVS = suppfree.build(base, definer=definer, ORDER=ORDER, FREE=FREE,
                                  modp=None)
U = set()
for c in BAD:
    m = suppfree.atom_supp(c, base, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
t0 = time.time()
cols = {}
for u in U:
    dv = {u: 1}
    for t in ORDER:
        a = definer[t]
        d = ad.dpart(a, t, bm)
        if d % P == 0:
            dv[t] = 0
            continue
        s = 0
        for w in L.avars[a]:
            if w == t:
                continue
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(a, w, bm) * dw
        dv[t] = (-s % P) * pow(d, -1, P) % P
    out = {}
    for c in CHECKS:
        s = 0
        for w in L.avars[c]:
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(c, w, bm) * dw
        if s % P:
            out[c] = s % P
    cols[u] = out
U = [u for u in U if cols[u]]
print(f'{len(U)} non-handle columns ({time.time()-t0:.0f}s)', flush=True)

atoms_moved = sorted(set().union(*[set(cols[u]) for u in U]))
EQS = sorted(set().union(*[set(L.atom2eq[a]) for a in atoms_moved]))
HOLD = [e for e in EQS if e not in set(FAILEQ)]
print(f'{len(atoms_moved)} atoms move -> {len(EQS)} equations touched, '
      f'{len(HOLD)} hold', flush=True)


def eqrow(e, u):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        d = cols[u].get(a, 0)
        if d:
            s += c * d
    return s % P


m = len(U)
A = [r for r in ([eqrow(e, u) for u in U] for e in HOLD) if any(r)]
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
print(f'rank {r_};  EQUATION-LEVEL KERNEL DIMENSION {len(KER)}', flush=True)

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
        fwd(v)
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
                    fwd(v)
                    s = L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))
                    if s > best:
                        best = s
                        T.save(v, os.path.join(HERE, 'EQK2_%d_%s.json' % (s, name)))
                out[label + 'score'] = best
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(DIRS)
sweep(tag, DIRS, evaluate, start, min(end, len(DIRS)), keyfn=lambda s: s[0],
      budget=540)
rs = load(tag)
if rs:
    print('\nbest directions at the deliverable:')
    for r in sorted(rs, key=lambda r: -max(r.get('safescore') or 0,
                                           r.get('fixscore') or 0))[:15]:
        print('   %-10s vary %-4d Gdeg %-4s safe-roots %-4s fix-roots %-4s '
              'score %s/%s' % (r['name'], r['vary'], r.get('Gdeg'),
                               r.get('saferoots'), r.get('fixroots'),
                               r.get('safescore'), r.get('fixscore')))
