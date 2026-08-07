"""S11 step 100: steer the two congruence targets to the cubic's solution.

close.py sets w3 and w2 to the cubic's root and A and B become EXACTLY ZERO -- the
A and B both vanish for the first time.  What it costs is the two advice
congruences whose left sides we just overwrote:

    a26731   x16742 = x19083        we set x16742 = w2*, so x19083 must follow
    a29539   x14853 = x1308         we set x14853 = w3*, so x1308  must follow

and valjac says those targets are not fixed: x19083 is moved by 170 free inputs and
x1308 by 79.  So the endgame is to steer them.

The complication is that the same inputs also move w1 = x12186 (179 of them), and the
cubic's root depends on w1 -- so this is a fixed point, not a one-shot solve.  Iterate:
solve the exact linear system for the current targets, apply, recompute the cubic from
the new w1 and w4, repeat.  Every step is measured, never predicted.

Usage: steer.py [state.json] [ROUNDS]
"""
import os, sys, json, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
import fpoly as F
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PIN_39013.json'
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 6
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
mv = json.load(open(os.path.join(HERE, 'valmovers.json')))
KN = sorted(set(mv['y1_tgt']) | set(mv['x2_tgt']))
KN = [u for u in KN if u not in (16742, 14853, 22162, 30213)]
print('%s: %d knobs move x19083 or x1308' % (src, len(KN)), flush=True)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-36s score %-6d A=%s B=%s  checks %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, nz), flush=True)
    return s, nz


def cubic_target(v):
    w1, w4 = v[12186] % P, v[24908] % P
    w5, w6, K = v[22162] % P, v[30213] % P, v[24453] % P
    S = (w4 + w6) % P
    cub = F.psub(F.pmul([(w5 + w1 + K) % P, 1],
                        [w5 * w5 % P, (-2 * w5) % P, 1]), [S * S % P])
    out = []
    for m in F.roots(cub):
        if (m - w5) % P == 0:
            continue
        w = S * ((m - w1) % P) % P * pow((m - w5) % P, -1, P) % P
        out.append((m, (w4 - w) % P))
    return out


def setcoords(v, m, Y1):
    v = list(v)
    v[14853] = (v[14853] // P) * P + m
    v[16742] = (v[16742] // P) * P + Y1
    ad.fwd(v, rounds=6)
    return v


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


def steer(v, tgt19083, tgt1308, knobs):
    """One exact linear step driving (x19083, x1308) to their targets."""
    b = [(tgt19083 - v[19083]) % P, (tgt1308 - v[1308]) % P]
    if not any(b):
        return v, True
    cols = []
    for u in knobs:
        w = list(v)
        w[u] = w[u] + 1
        ad.fwd(w, rounds=6)
        cols.append([(w[19083] - v[19083]) % P, (w[1308] - v[1308]) % P])
    n, m = 2, len(knobs)
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
    ok = (w[19083] % P == tgt19083 % P and w[1308] % P == tgt1308 % P)
    return w, ok


best, bestv = report(v, 'start')[0], list(v)
t0 = time.time()
for rnd in range(ROUNDS):
    tg = cubic_target(v)
    if not tg:
        print('round %d: cubic has no usable root' % rnd, flush=True)
        break
    m, Y1 = tg[0]
    print('\nround %d: cubic root w3* = %s..., w2* = %s...'
          % (rnd, str(m)[:22], str(Y1)[:22]), flush=True)
    for K in (12, 30, len(KN)):
        w, ok = steer(v, Y1, m, KN[:K])
        if not ok:
            print('   %d knobs: could not hit both targets' % K, flush=True)
            continue
        w = setcoords(w, m, Y1)
        s, nz = report(w, '   %d knobs: steered + set' % K)
        w = lift(w)
        s, nz = report(w, '   %d knobs: after lift' % K)
        if s > best:
            best, bestv = s, list(w)
            T.save(w, os.path.join(HERE, 'ST_%d.json' % s))
            print('      *** saved ST_%d.json (%.0fs)' % (s, time.time() - t0),
                  flush=True)
        if s >= best:
            v = w
    if not tg:
        break
print('\nbest %d' % best)
T.save(bestv, os.path.join(HERE, 'ST_best_%d.json' % best))
