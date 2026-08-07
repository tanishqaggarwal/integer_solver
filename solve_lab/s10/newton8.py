"""S11 step 102: the coordinate map has RANK 8 -- solve all eight at once, by Newton.

coordjac's data answers the question that matters.  Over the 264 free inputs with any
effect, the map to

    (x1, y1, x2, y2, x3, y3, x19083, x1308)   plus A and B

has **rank 8 of 10**, and the only two relations are the linearisations of A and B
themselves:

    0 = c1*x1 + c2*y1 + c3*x2 + c4*y2 + c5*x3 + A
    0 = d1*x1 + d2*y1 + d3*x2 + d4*y2 + d5*x3 + d6*y3 + B

So the eight coordinates are INDEPENDENTLY STEERABLE to first order.  close2 failed
only because it used the 180 columns that move x1 or x19083; the full 264 has the
rank.

Target: keep x3 = C2 and y3 = C1' (so a1618 and a688 hold), keep x2 and y2, and move
x1 and y1 to the pair-(x1,y1) LINEAR solution of A = B = 0, dragging x19083 with y1 so
a26731 survives and leaving x1308 alone so a29539 does.  Eight targets, 264 knobs,
solved with the exact first-order columns and then Newton-iterated because the map is
a genuine polynomial.

Usage: newton8.py [state.json] [ROUNDS]
"""
import os, sys, json, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PIN_39013.json'
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 8
KEY = [12186, 16742, 14853, 24908, 22162, 30213, 19083, 1308]
NM = dict(zip(KEY, ['x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'y1t', 'x2t']))
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
rs = [json.loads(l) for l in
      open(os.path.join(HERE, 'runs', 'coord_PIN_39013.jsonl'))
      if l.strip().endswith('}')]
rs = [r.get('result', r) for r in rs]
KN = sorted({r['u'] for r in rs if r['d']})
print('%s: %d free inputs with any effect on the coordinates' % (src, len(KN)),
      flush=True)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-34s score %-6d A=%-5s B=%-5s checks %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, nz), flush=True)
    return s, nz


def targets(v):
    """Keep x2, y2, x3, y3; move x1, y1 to the linear solution; drag x19083."""
    x2, y2 = v[14853] % P, v[24908] % P
    x3, y3, K = v[22162] % P, v[30213] % P, v[24453] % P
    d = (x2 - x3) % P
    x1s = (((y2 + y3) % P) ** 2 % P * pow(d * d % P, -1, P) - x3 - x2 - K) % P
    y1s = (((x1s - x3) % P) * y2 - y3 * ((x2 - x1s) % P)) % P * pow(d, -1, P) % P
    return {12186: x1s, 16742: y1s, 19083: y1s, 14853: x2, 24908: y2,
            22162: x3, 30213: y3, 1308: v[1308] % P}


def jac(v, knobs, keys):
    vm = [x % P for x in v]
    cols = []
    for u in knobs:
        dv = {u: 1}
        for t in ad.ORDER:
            a = L.definer[t]
            dd = ad.dpart(a, t, vm)
            if dd % P == 0:
                dv[t] = 0
                continue
            s = 0
            for w in L.avars[a]:
                if w == t:
                    continue
                dw = dv.get(w, 0)
                if dw:
                    s += ad.dpart(a, w, vm) * dw
            dv[t] = (-s % P) * pow(dd, -1, P) % P
        cols.append([dv.get(k, 0) % P for k in keys])
    return cols


def step(v, tg, knobs):
    keys = list(tg)
    b = [(tg[k] - v[k]) % P for k in keys]
    if not any(b):
        return v, True
    cols = jac(v, knobs, keys)
    n, m = len(keys), len(knobs)
    A = [[cols[j][i] for j in range(m)] + [b[i]] for i in range(n)]
    piv, r_ = [], 0
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
    if any(A[i][m] for i in range(r_, n)):
        return v, False
    d = [0] * m
    for i, j in enumerate(piv):
        d[j] = A[i][m]
    w = list(v)
    for j, u in enumerate(knobs):
        if d[j]:
            w[u] = w[u] + d[j]
    ad.fwd(w, rounds=6)
    return w, True


def lift(v, budget=25):
    _, fl, SV = suppfree.build(v, modp=None)
    for _ in range(budget):
        av = L.all_atom_values(v)
        todo = [a for a in range(L.NA) if a not in L.atom_out and av[a]
                and av[a] % P == 0]
        cur = L.NEQ - len(L.failing_eqs(av))
        moved = False
        for c in todo:
            mm = suppfree.atom_supp(c, v, SV, modp=None)
            for i in range(len(fl)):
                if not ((mm >> i) & 1):
                    continue
                t = fl[i]
                g = jacZ(t, v, [c]).get(c, 0)
                if not g or g % P or av[c] % g:
                    continue
                w = list(v)
                w[t] = w[t] - av[c] // g
                ad.fwd(w, rounds=6)
                a2 = L.all_atom_values(w)
                if a2[c] == 0 and L.NEQ - len(L.failing_eqs(a2)) >= cur:
                    v, moved = w, True
                    break
            if moved:
                break
        if not moved:
            break
    return v


best, bestv = report(v, 'start')[0], list(v)
tg = targets(v)
print('targets: ' + ', '.join('%s=%s...' % (NM[k], str(tg[k])[:14]) for k in tg),
      flush=True)
t0 = time.time()
cur = list(v)
for rnd in range(ROUNDS):
    tg = targets(cur)
    w, ok = step(cur, tg, KN)
    if not ok:
        print('round %d: system inconsistent over %d knobs' % (rnd, len(KN)),
              flush=True)
        break
    err = max(abs((w[k] - tg[k]) % P) for k in tg)
    s, nz = report(w, 'round %d: stepped' % rnd)
    print('     residual on the 8 targets: %s'
          % ('ZERO' if err == 0 else 'nonzero'), flush=True)
    w2 = lift(w)
    s, nz = report(w2, 'round %d: lifted' % rnd)
    if s > best:
        best, bestv = s, list(w2)
        T.save(w2, os.path.join(HERE, 'N8_%d.json' % s))
        print('     *** saved N8_%d.json (%.0fs)' % (s, time.time() - t0),
              flush=True)
    cur = w
print('\nbest %d' % best)
T.save(bestv, os.path.join(HERE, 'N8_best_%d.json' % best))
