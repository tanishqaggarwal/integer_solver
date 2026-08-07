"""S11 step 104: the w-map is DIAGONAL -- one advice value per w.

Of the 264 free inputs that move the six values, only EIGHT are non-boolean, and
they form a diagonal system -- each w has exactly one non-boolean knob:

    w1  = x12186 <- x22649       w2  = x16742 <- itself
    w3  = x14853 <- itself       w4  = x24908 <- x31339
    w5  = x22162 <- itself       w6  = x30213 <- itself
    w2t = x19083 <- x8778        w3t = x1308  <- x6418

(the other 256 movers are boolean, and driving a boolean to an arbitrary residue
breaks its own b^2 = b constraint, which is why the 264-knob Newton cost 91.)

Each of those eight knobs carries exactly one congruence, so freeing one w
costs exactly one congruence:

    w1  free  <=>  break a2423  (x22649 = x29524)      w2 free <=> break a26731
    w4  free  <=>  break a33796 (x31339 = x6858)       w3 free <=> break a29539
    w2t free  <=>  break a33929 (x8778  = x16144)      w5, w6 <=> break a1618/a688

and A = B = 0 needs exactly TWO values.  Two of the six pairs solve LINEARLY --
(w1,w2) and (w3,w4) -- one needs a cubic, and (w5,w6) solves linearly.  So price every pair:
free its two values, set them to the exact solution, lift, and measure.

Usage: pairfix.py [state.json]
"""
import os, sys, time
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
v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v0, rounds=6)
KNOB = {12186: 22649, 16742: 16742, 14853: 14853, 24908: 31339,
        22162: 22162, 30213: 30213}
NM = {12186: 'w1', 16742: 'w2', 14853: 'w3', 24908: 'w4',
      22162: 'w5', 30213: 'w6'}


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-40s score %-6d A=%-5s B=%-5s checks %d %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, len(nz), nz[:12]),
          flush=True)
    return s, nz


def setcoord(v, coord, target):
    """Drive one value to a residue via its single non-boolean knob."""
    u = KNOB[coord]
    if u == coord:
        v = list(v)
        v[coord] = (v[coord] // P) * P + target % P
        ad.fwd(v, rounds=6)
        return v
    w = list(v)
    w[u] = w[u] + 1
    ad.fwd(w, rounds=6)
    g = (w[coord] - v[coord]) % P
    if g == 0:
        return None
    d = (target - v[coord]) % P * pow(g, -1, P) % P
    w = list(v)
    w[u] = w[u] + d
    ad.fwd(w, rounds=6)
    return w if w[coord] % P == target % P else None


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


def coords(v):
    return (v[12186] % P, v[16742] % P, v[14853] % P, v[24908] % P,
            v[22162] % P, v[30213] % P, v[24453] % P)


def sols(pair, c):
    """Exact solutions of A = B = 0 for a pair, others held."""
    w1, w2, w3, w4, w5, w6, K = c
    out = []
    if pair == ('w1', 'w2'):
        d = (w3 - w5) % P
        if d:
            w1s = (((w4 + w6) % P) ** 2 % P * pow(d * d % P, -1, P)
                   - w5 - w3 - K) % P
            w2s = (((w1s - w5) % P) * w4 - w6 * ((w3 - w1s) % P)) % P \
                * pow(d, -1, P) % P
            out.append({12186: w1s, 16742: w2s})
    elif pair == ('w3', 'w4'):
        d = (w1 - w5) % P
        if d:
            w3s = (((w6 + w2) % P) ** 2 % P * pow(d * d % P, -1, P)
                   - w5 - w1 - K) % P
            w4s = (w2 + ((w6 + w2) % P) * ((w3s - w1) % P) % P
                   * pow(d, -1, P)) % P
            out.append({14853: w3s, 24908: w4s})
    elif pair == ('w5', 'w6'):
        d = (w3 - w1) % P
        if d:
            w5s = (((w4 - w2) % P) ** 2 % P * pow(d * d % P, -1, P)
                   - w1 - w3 - K) % P
            w6s = (((w1 - w5s) % P) * ((w4 - w2) % P) % P
                   * pow(d, -1, P) - w2) % P
            out.append({22162: w5s, 30213: w6s})
    elif pair == ('w3', 'w2'):
        S = (w4 + w6) % P
        cub = F.psub(F.pmul([(w5 + w1 + K) % P, 1],
                            [w5 * w5 % P, (-2 * w5) % P, 1]), [S * S % P])
        for m in F.roots(cub):
            if (m - w5) % P == 0:
                continue
            w = S * ((m - w1) % P) % P * pow((m - w5) % P, -1, P) % P
            out.append({14853: m, 16742: (w4 - w) % P})
    elif pair == ('w1', 'w4'):
        S = (w6 + w2) % P
        cub = F.psub(F.pmul([(w5 + w3 + K) % P, 1],
                            [w5 * w5 % P, (-2 * w5) % P, 1]), [S * S % P])
        for m in F.roots(cub):
            if (m - w5) % P == 0:
                continue
            w = S * ((w3 - m) % P) % P * pow((m - w5) % P, -1, P) % P
            out.append({12186: m, 24908: (w2 + w) % P})
    return out


base = report(v0, 'start (%s)' % src)[0]
best, bestv = base, list(v0)
t0 = time.time()
for pair in [('w1', 'w2'), ('w3', 'w4'), ('w5', 'w6'), ('w3', 'w2'), ('w1', 'w4')]:
    c = coords(v0)
    ss = sols(pair, c)
    print('\n=== pair %s : %d exact solution(s) ===' % (str(pair), len(ss)),
          flush=True)
    for sol in ss:
        v = list(v0)
        ok = True
        for coord, tgt in sol.items():
            w = setcoord(v, coord, tgt)
            if w is None:
                print('   could not drive %s' % NM[coord], flush=True)
                ok = False
                break
            v = w
        if not ok:
            continue
        s, nz = report(v, '   set %s' % '+'.join(NM[k] for k in sol))
        v = lift(v)
        s, nz = report(v, '   lifted')
        T.save(v, os.path.join(HERE, 'PF_%s_%d.json'
                               % ('_'.join(pair), s)))
        if s > best:
            best, bestv = s, list(v)
            T.save(v, os.path.join(HERE, 'PF_best_%d.json' % s))
            print('      *** NEW BEST %d  (%.0fs)' % (s, time.time() - t0),
                  flush=True)
print('\nbest %d' % best)
