"""S11 step 103: Newton on the coordinates, then the advice DAG, then the lift -- looped.

newton8 converges: by round 1 the residual on all eight coordinate targets is ZERO,
A and B are both exactly zero, and a26731/a29539 hold because x19083 and x1308 were
dragged along as targets.  **The two-condition primitive closes with its congruences intact.**

What it costs is the OTHER advice congruences -- a33929 (x8778 = x16144), a2423,
a21617, a7930, a33796 and their bundles -- because the 264 knobs feed those chains
too.  But those are exactly what advgraph.py solves in one topological sweep, since
each of those advice values is free and only has to follow its target.

So compose the three exact operations and iterate:

    Newton-step the eight coordinates  ->  Gauss-Seidel the advice DAG  ->  integer lift

Each is exact and measured; the loop is the only thing that is heuristic.

Usage: full.py [state.json] [ROUNDS]
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
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
rs = [json.loads(l) for l in
      open(os.path.join(HERE, 'runs', 'coord_PIN_39013.jsonl'))
      if l.strip().endswith('}')]
KN = sorted({r.get('result', r)['u'] for r in rs if r.get('result', r)['d']})
TWO = [(8778, 16144), (14623, 27522), (16742, 19083), (22649, 29524),
       (24548, 25442), (31339, 6858), (14853, 1308)]
ADVORDER = [8778, 24548, 14623, 16742, 22649, 31339, 14853]


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-30s score %-6d A=%-5s B=%-5s checks %d %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, len(nz), nz[:14]),
          flush=True)
    return s, nz


def targets(v):
    x2, y2 = v[14853] % P, v[24908] % P
    x3, y3, K = v[22162] % P, v[30213] % P, v[24453] % P
    d = (x2 - x3) % P
    x1s = (((y2 + y3) % P) ** 2 % P * pow(d * d % P, -1, P) - x3 - x2 - K) % P
    y1s = (((x1s - x3) % P) * y2 - y3 * ((x2 - x1s) % P)) % P * pow(d, -1, P) % P
    return {12186: x1s, 16742: y1s, 19083: y1s, 14853: x2, 24908: y2,
            22162: x3, 30213: y3, 1308: v[1308] % P}


def step(v, tg, knobs):
    keys = list(tg)
    b = [(tg[k] - v[k]) % P for k in keys]
    if not any(b):
        return v, True
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


def advice(v, hold=()):
    v = list(v)
    for _ in range(3):
        for t in ADVORDER:
            if t in hold:
                continue
            y = dict(TWO)[t]
            v[t] = (v[t] // P) * P + v[y] % P
            ad.fwd(v, rounds=6)
    return v


def lift(v, budget=30):
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
t0 = time.time()
cur = list(v)
for rnd in range(ROUNDS):
    w, ok = step(cur, targets(cur), KN)
    if not ok:
        print('round %d: inconsistent' % rnd, flush=True)
        break
    report(w, 'round %d: newton' % rnd)
    w = advice(w, hold=(16742, 14853))
    report(w, 'round %d: advice' % rnd)
    w = lift(w)
    s, nz = report(w, 'round %d: lift' % rnd)
    if s > best:
        best, bestv = s, list(w)
        T.save(w, os.path.join(HERE, 'FULL_%d.json' % s))
        print('   *** saved FULL_%d.json (%.0fs)' % (s, time.time() - t0), flush=True)
    T.save(w, os.path.join(HERE, 'FULL_r%d_%d.json' % (rnd, s)))
    cur = w
print('\nbest %d' % best)
T.save(bestv, os.path.join(HERE, 'FULL_best_%d.json' % best))
