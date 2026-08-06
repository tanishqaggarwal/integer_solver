"""Augmented model: single-variable knobs PLUS compound directions (integer combinations
of variable moves whose disturbance of the atoms outside A cancels exactly).
Then the exact integer optimisation over S13 and realisation of the winner."""
import pickle, collections, sys, time, itertools, json
import lib as L, model as MD, opt, ihs
from fractions import Fraction

v0 = opt.init()
P = L.P
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])
A = set(MD.confined_atoms(S13))
base = MD.build(S13, v0, verbose=False)
C2 = MD.load_census('24')
rev = C2['rev']

MAXEXTRA = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 2

E = set()
for _ in range(ROUNDS):
    pool = set()
    for a in A | E:
        pool |= rev.get(a, set())
    newE = set()
    for x in sorted(pool):
        v1, t1 = MD.move(v0, {x: v0[x] + 1}, A)
        ex = set(t1) - A - E
        if 0 < len(ex) <= MAXEXTRA:
            newE |= ex
    if newE <= E:
        break
    E |= newE
B = A | E
pool = set()
for a in B:
    pool |= rev.get(a, set())
vars_ = []
cols = {}
for x in sorted(pool):
    v1, t1 = MD.move(v0, {x: v0[x] + 1}, A)
    if not t1 or not set(t1) <= B:
        continue
    v5, t5 = MD.move(v0, {x: v0[x] + 5}, A)
    if set(t5) != set(t1) or any(t5[a] - MD.BASEP[a] != 5 * (t1[a] - MD.BASEP[a]) for a in t1):
        continue
    vars_.append(x)
    cols[x] = {a: t1[a] - MD.BASEP[a] for a in t1}
print(f'|E|={len(E)}  variables usable = {len(vars_)}')

Eidx = sorted(E)
nc = len(vars_)
m = [[Fraction(cols[x].get(a, 0)) for x in vars_] for a in Eidx]
piv = []
r = 0
for c in range(nc):
    p = next((i for i in range(r, len(m)) if m[i][c]), None)
    if p is None:
        continue
    m[r], m[p] = m[p], m[r]
    pv = m[r][c]
    for i in range(len(m)):
        if i != r and m[i][c]:
            f = m[i][c] / pv
            for j in range(c, nc):
                m[i][j] -= f * m[r][j]
    piv.append(c)
    r += 1
free = [c for c in range(nc) if c not in piv]
import math
basis = []
for fc in free:
    vec = [Fraction(0)] * nc
    vec[fc] = Fraction(1)
    for i, pc in enumerate(piv):
        vec[pc] = -m[i][fc] / m[i][pc]
    den = 1
    for q in vec:
        den = den * q.denominator // math.gcd(den, q.denominator)
    iv = [int(q * den) for q in vec]
    g = 0
    for q in iv:
        g = math.gcd(g, abs(q))
    if g > 1:
        iv = [q // g for q in iv]
    basis.append(iv)
print(f'kernel dim = {len(basis)}')

# --- augmented knob list: name -> (variable-move dict, A-motion dict) -----------
knobs = []
for x in base['knobs']:
    knobs.append((('v', x), {x: 1}, {a: base['M'][a].get(x, 0) for a in A if base['M'][a].get(x, 0)}))
for bi, b in enumerate(basis):
    mv = {vars_[j]: b[j] for j in range(nc) if b[j]}
    mo = collections.defaultdict(int)
    for x, c in mv.items():
        for a, d in cols[x].items():
            mo[a] += c * d
    mo = {a: c for a, c in mo.items() if c and a in A}
    if not mo:
        continue
    if any(a not in A for a in mo):
        continue
    knobs.append((('c', bi), mv, mo))
print(f'augmented knobs = {len(knobs)}  (singles {len(base["knobs"])}, compounds {len(knobs)-len(base["knobs"])})')

names = [k[0] for k in knobs]
mod = {'knobs': names, 'A': sorted(A), 'S': sorted(S13),
       'M': {a: {names[i]: knobs[i][2].get(a, 0) for i in range(len(knobs)) if knobs[i][2].get(a, 0)} for a in A}}
conds = []
for i in sorted(S13):
    mm, sq, co = L.eq_atoms[i]
    lin = [(a, c) for a, c in co.items() if not MD.IS_SQ[a]]
    sqs = [(a, c) for a, c in co.items() if MD.IS_SQ[a]]
    rows = []

    def linrow(terms):
        row = {n: 0 for n in names}
        rhs = 0
        for a, c in terms:
            rhs -= c * MD.BASEP[a]
            for n, mc in mod['M'].get(a, {}).items():
                row[n] += c * mc
        return row, rhs
    if not sqs:
        rows.append(linrow(lin))
    elif not lin and len(sqs) == 1:
        a, c = sqs[0]
        rows.append(({n: mod['M'].get(a, {}).get(n, 0) for n in names}, -MD.BASEP[a]))
    else:
        for a, c in sqs:
            rows.append(({n: mod['M'].get(a, {}).get(n, 0) for n in names}, -MD.BASEP[a]))
        rows.append(linrow(lin))
    conds.append((i, rows))
mod['conds'] = conds
print('D =', opt.rank_of(mod))

t0 = time.time()
sz, T, z, calls = opt.maximise(mod, S13)
print(f'MAX zeroable with compound knobs = {sz}  {T}  -> failing {13-sz}   ({time.time()-t0:.0f}s)')

if sz > 6:
    seeds = collections.defaultdict(int)
    for i, n in enumerate(names):
        d = z.get(n, 0)
        if d:
            for x, c in knobs[i][1].items():
                seeds[x] += d * c
    sd = {x: v0[x] + c for x, c in seeds.items() if c}
    v, tou = MD.move(v0, sd, A)
    av = L.all_atom_values(v)
    fails = L.failing_eqs(av)
    print('REALISED failing:', len(fails), fails)
    if len(fails) <= 6:
        L.save(v, '/home/user/integer_solver/solve_lab/s9/eff/cand6.json')
        print('saved cand6.json')
