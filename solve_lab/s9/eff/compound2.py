"""General compound-knob search: allow a set E of extra atoms to be disturbed but force
them back to zero by a simultaneous integer combination of variable moves.  Any resulting
direction outside the current knob lattice would break one of the two surviving congruences.
"""
import pickle, collections, sys, time
import lib as L, model as MD, opt
from fractions import Fraction

v0 = opt.init()
P = L.P
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])
A = set(MD.confined_atoms(S13))
mod = MD.build(S13, v0, verbose=False)
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
print(f'extra atom pool E: {len(E)}')

B = A | E
pool = set()
for a in B:
    pool |= rev.get(a, set())
knobs = []
cols = {}
for x in sorted(pool):
    v1, t1 = MD.move(v0, {x: v0[x] + 1}, A)
    if not t1 or not set(t1) <= B:
        continue
    v5, t5 = MD.move(v0, {x: v0[x] + 5}, A)
    if set(t5) != set(t1) or any(t5[a] - MD.BASEP[a] != 5 * (t1[a] - MD.BASEP[a]) for a in t1):
        continue
    knobs.append(x)
    cols[x] = {a: t1[a] - MD.BASEP[a] for a in t1}
print(f'variables with footprint inside A u E: {len(knobs)}')

Eidx = sorted(E)
# rational kernel of the |E| x |knobs| matrix, then look at the induced motion on A
rows = [[Fraction(cols[x].get(a, 0)) for x in knobs] for a in Eidx]
m = [r[:] for r in rows]
nr, nc = len(m), len(knobs)
piv = []
r = 0
for c in range(nc):
    p = next((i for i in range(r, nr) if m[i][c]), None)
    if p is None:
        continue
    m[r], m[p] = m[p], m[r]
    pv = m[r][c]
    for i in range(nr):
        if i != r and m[i][c]:
            f = m[i][c] / pv
            for j in range(c, nc):
                m[i][j] -= f * m[r][j]
    piv.append(c)
    r += 1
free = [c for c in range(nc) if c not in piv]
print(f'rank(E-block) = {r}, kernel dim = {len(free)}')

# basis of the kernel
basis = []
for fc in free:
    vec = [Fraction(0)] * nc
    vec[fc] = Fraction(1)
    for i, pc in enumerate(piv):
        vec[pc] = -m[i][fc] / m[i][pc]
    den = 1
    for q in vec:
        den = den * q.denominator // __import__('math').gcd(den, q.denominator)
    basis.append([int(q * den) for q in vec])

Aidx = sorted(a for a in A if a != 37887)
rowsA = [[cols[x].get(a, 0) for x in knobs] for a in Aidx]
cur = [[mod['M'][a].get(x, 0) for x in mod['knobs']] for a in Aidx]


def snf_diag(Ain):
    Am = [r[:] for r in Ain]
    mm = len(Am)
    nn = len(Am[0]) if Am and Am[0] else 0
    d = []
    t = 0
    while t < min(mm, nn):
        pv = None
        for i in range(t, mm):
            for j in range(t, nn):
                if Am[i][j]:
                    pv = (i, j)
                    break
            if pv:
                break
        if pv is None:
            break
        i, j = pv
        Am[t], Am[i] = Am[i], Am[t]
        for rr in Am:
            rr[t], rr[j] = rr[j], rr[t]
        while True:
            ch = False
            for i in range(t + 1, mm):
                if Am[i][t]:
                    q = Am[i][t] // Am[t][t]
                    for k in range(nn):
                        Am[i][k] -= q * Am[t][k]
                    if Am[i][t]:
                        Am[t], Am[i] = Am[i], Am[t]
                        ch = True
            for j in range(t + 1, nn):
                if Am[t][j]:
                    q = Am[t][j] // Am[t][t]
                    for rr in Am:
                        rr[j] -= q * rr[t]
                    if Am[t][j]:
                        for rr in Am:
                            rr[t], rr[j] = rr[j], rr[t]
                        ch = True
            if ch:
                continue
            if all(Am[i][t] == 0 for i in range(t + 1, mm)) and all(Am[t][j] == 0 for j in range(t + 1, nn)):
                break
        d.append(abs(Am[t][t]))
        t += 1
    return d


gen = [list(c) for c in zip(*cur)]                       # current knob columns on A
newdirs = 0
for b in basis:
    col = [sum(rowsA[i][j] * b[j] for j in range(nc)) for i in range(len(Aidx))]
    if any(col):
        gen.append(col)
        newdirs += 1
print(f'compound directions produced: {newdirs}')
Mall = [list(c) for c in zip(*gen)]
d_old = snf_diag(cur)
d_new = snf_diag(Mall)
print('current lattice invariant factors :', ['P' if x == P else ('7376877P' if x == 7376877 * P else x) for x in d_old])
print('with compound knobs               :', ['P' if x == P else ('7376877P' if x == 7376877 * P else x) for x in d_new])
print('C before =', sum(1 for x in d_old if x and x % P == 0),
      '  C after =', sum(1 for x in d_new if x and x % P == 0))
