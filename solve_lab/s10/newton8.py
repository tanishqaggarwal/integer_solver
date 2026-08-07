"""S11 step 102: the w-map has RANK 8 -- solve all eight at once, by Newton.

valjac's data answers the question that matters.  Over the 264 free inputs with any
effect, the map to

    (w1, w2, w3, w4, w5, w6, x19083, x1308)   plus A and B

has **rank 8 of 10**, and the only two relations are the linearisations of A and B
themselves:

    0 = c1*w1 + c2*w2 + c3*w3 + c4*w4 + c5*w5 + A
    0 = d1*w1 + d2*w2 + d3*w3 + d4*w4 + d5*w5 + d6*w6 + B

So the eight values are INDEPENDENTLY STEERABLE to first order.  close2 failed
only because it used the 180 columns that move w1 or x19083; the full 264 has the
rank.

Target: keep w5 = C2 and w6 = C1' (so a1618 and a688 hold), keep w3 and w4, and move
w1 and w2 to the pair-(w1,w2) LINEAR solution of A = B = 0, dragging x19083 with w2 so
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
NM = dict(zip(KEY, ['w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w2t', 'w3t']))
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
rs = [json.loads(l) for l in
      open(os.path.join(HERE, 'runs', 'coord_PIN_39013.jsonl'))
      if l.strip().endswith('}')]
rs = [r.get('result', r) for r in rs]
KN = sorted({r['u'] for r in rs if r['d']})
print('%s: %d free inputs with any effect on the values' % (src, len(KN)),
      flush=True)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-34s score %-6d A=%-5s B=%-5s checks %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, nz), flush=True)
    return s, nz


def targets(v):
    """Keep w3, w4, w5, w6; move w1, w2 to the linear solution; drag x19083."""
    w3, w4 = v[14853] % P, v[24908] % P
    w5, w6, K = v[22162] % P, v[30213] % P, v[24453] % P
    d = (w3 - w5) % P
    w1s = (((w4 + w6) % P) ** 2 % P * pow(d * d % P, -1, P) - w5 - w3 - K) % P
    w2s = (((w1s - w5) % P) * w4 - w6 * ((w3 - w1s) % P)) % P * pow(d, -1, P) % P
    return {12186: w1s, 16742: w2s, 19083: w2s, 14853: w3, 24908: w4,
            22162: w5, 30213: w6, 1308: v[1308] % P}


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
