"""S11 step 101: solve the addition in the pair (x1, y1) -- which is LINEAR -- and
steer x19083 with it so the congruence survives.

Eliminating y1 between A and B gives, remarkably,

    B  =>  y1 = [(x1-x3)*y2 - y3*(x2-x1)] / (x2-x3)
    A  =>  (x3 + x1 + x2 + K) * (x2-x3)^2 = (y2+y3)^2

and the second is **LINEAR in x1**.  So the pair (x1, y1) always has a solution, with
no cubic and no root-existence question:

    x1* = (y2+y3)^2/(x2-x3)^2 - x3 - x2 - K            y1* from B

Two of the four pairs behave this way -- (x1,y1) and (x2,y2) -- while (x2,y1) needs a
cubic and (x3,y3) is the pair already used.

The reason this pair is the right one: y1 = x16742 is FREE, so y1* can just be
written in, and its congruence a26731 (x16742 = x19083) is repaired in the same move,
because x19083 and x1 = x12186 are moved by almost the same free inputs (170 and 179,
overlapping).  Steer BOTH at once -- x12186 to x1* and x19083 to y1* -- and A = B = 0
while a26731 still holds and a29539 is never touched.

Usage: close2.py [state.json] [ROUNDS] [NKNOB]
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
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
NK = int(sys.argv[3]) if len(sys.argv) > 3 else 60
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
mv = json.load(open(os.path.join(HERE, 'coordmovers.json')))
KN = [u for u in sorted(set(mv['x1']) | set(mv['y1_tgt']))
      if u not in (16742, 14853, 22162, 30213)]
print('%s: %d knobs move x12186 or x19083' % (src, len(KN)), flush=True)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-38s score %-6d A=%-5s B=%-5s checks %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, nz), flush=True)
    return s, nz


def target(v):
    """x1*, y1* from the LINEAR elimination."""
    x2, y2 = v[14853] % P, v[24908] % P
    x3, y3, K = v[22162] % P, v[30213] % P, v[24453] % P
    d = (x2 - x3) % P
    if d == 0:
        return None
    x1s = ((y2 + y3) % P) ** 2 % P * pow(d * d % P, -1, P) % P
    x1s = (x1s - x3 - x2 - K) % P
    y1s = (((x1s - x3) % P) * y2 - y3 * ((x2 - x1s) % P)) % P * pow(d, -1, P) % P
    return x1s, y1s


def multi_steer(v, targets, knobs):
    """Exact linear step driving several variables to given residues."""
    b = [(t - v[u]) % P for u, t in targets]
    if not any(b):
        return v, True
    cols = []
    for k in knobs:
        w = list(v)
        w[k] = w[k] + 1
        ad.fwd(w, rounds=6)
        cols.append([(w[u] - v[u]) % P for u, _ in targets])
    n, m = len(targets), len(knobs)
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
    for j, k in enumerate(knobs):
        if d[j]:
            w[k] = w[k] + d[j]
    ad.fwd(w, rounds=6)
    return w, all(w[u] % P == t % P for u, t in targets)


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
t0 = time.time()
cur = list(v)
for rnd in range(ROUNDS):
    tg = target(cur)
    if tg is None:
        print('round %d: x2 == x3, degenerate' % rnd)
        break
    x1s, y1s = tg
    print('\nround %d: x1* = %s...   y1* = %s...'
          % (rnd, str(x1s)[:24], str(y1s)[:24]), flush=True)
    for K in (8, 20, min(NK, len(KN)), len(KN)):
        w, ok = multi_steer(cur, [(12186, x1s), (19083, y1s)], KN[:K])
        if not ok:
            print('   %-4d knobs: targets not reachable' % K, flush=True)
            continue
        w[16742] = (w[16742] // P) * P + y1s % P
        ad.fwd(w, rounds=6)
        s, nz = report(w, '   %d knobs: steered' % K)
        w = lift(w)
        s, nz = report(w, '   %d knobs: lifted' % K)
        if s > best:
            best, bestv = s, list(w)
            T.save(w, os.path.join(HERE, 'C2_%d.json' % s))
            print('      *** saved C2_%d.json  (%.0fs)' % (s, time.time() - t0),
                  flush=True)
        if s >= report(cur, '   (current)')[0]:
            cur = w
print('\nbest %d' % best)
T.save(bestv, os.path.join(HERE, 'C2_best_%d.json' % best))
