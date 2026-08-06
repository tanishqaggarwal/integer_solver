"""S11 step 104: the coordinate map is DIAGONAL -- one advice value per coordinate.

Of the 264 free inputs that move the coordinates, only EIGHT are non-boolean, and
they form a diagonal system -- each coordinate has exactly one non-boolean knob:

    x1  = x12186 <- x22649       y1  = x16742 <- itself
    x2  = x14853 <- itself       y2  = x24908 <- x31339
    x3  = x22162 <- itself       y3  = x30213 <- itself
    y1t = x19083 <- x8778        x2t = x1308  <- x6418

(the other 256 movers are boolean, and driving a boolean to an arbitrary residue
breaks its own b^2 = b constraint, which is why the 264-knob Newton cost 91.)

Each of those eight knobs carries exactly one congruence, so freeing one coordinate
costs exactly one congruence:

    x1  free  <=>  break a2423  (x22649 = x29524)      y1 free <=> break a26731
    y2  free  <=>  break a33796 (x31339 = x6858)       x2 free <=> break a29539
    y1t free  <=>  break a33929 (x8778  = x16144)      x3, y3 <=> break a1618/a688

and A = B = 0 needs exactly TWO coordinates.  Two of the six pairs solve LINEARLY --
(x1,y1) and (x2,y2) -- one needs a cubic, and (x3,y3) is ecfix.  So price every pair:
free its two coordinates, set them to the exact solution, lift, and measure.

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
NM = {12186: 'x1', 16742: 'y1', 14853: 'x2', 24908: 'y2',
      22162: 'x3', 30213: 'y3'}


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-40s score %-6d A=%-5s B=%-5s checks %d %s'
          % (tag, s, v[35389] % P == 0, v[6671] % P == 0, len(nz), nz[:12]),
          flush=True)
    return s, nz


def setcoord(v, coord, target):
    """Drive one coordinate to a residue via its single non-boolean knob."""
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
    x1, y1, x2, y2, x3, y3, K = c
    out = []
    if pair == ('x1', 'y1'):
        d = (x2 - x3) % P
        if d:
            x1s = (((y2 + y3) % P) ** 2 % P * pow(d * d % P, -1, P)
                   - x3 - x2 - K) % P
            y1s = (((x1s - x3) % P) * y2 - y3 * ((x2 - x1s) % P)) % P \
                * pow(d, -1, P) % P
            out.append({12186: x1s, 16742: y1s})
    elif pair == ('x2', 'y2'):
        d = (x1 - x3) % P
        if d:
            x2s = (((y3 + y1) % P) ** 2 % P * pow(d * d % P, -1, P)
                   - x3 - x1 - K) % P
            y2s = (y1 + ((y3 + y1) % P) * ((x2s - x1) % P) % P
                   * pow(d, -1, P)) % P
            out.append({14853: x2s, 24908: y2s})
    elif pair == ('x3', 'y3'):
        d = (x2 - x1) % P
        if d:
            x3s = (((y2 - y1) % P) ** 2 % P * pow(d * d % P, -1, P)
                   - x1 - x2 - K) % P
            y3s = (((x1 - x3s) % P) * ((y2 - y1) % P) % P
                   * pow(d, -1, P) - y1) % P
            out.append({22162: x3s, 30213: y3s})
    elif pair == ('x2', 'y1'):
        S = (y2 + y3) % P
        cub = F.psub(F.pmul([(x3 + x1 + K) % P, 1],
                            [x3 * x3 % P, (-2 * x3) % P, 1]), [S * S % P])
        for m in F.roots(cub):
            if (m - x3) % P == 0:
                continue
            w = S * ((m - x1) % P) % P * pow((m - x3) % P, -1, P) % P
            out.append({14853: m, 16742: (y2 - w) % P})
    elif pair == ('x1', 'y2'):
        S = (y3 + y1) % P
        cub = F.psub(F.pmul([(x3 + x2 + K) % P, 1],
                            [x3 * x3 % P, (-2 * x3) % P, 1]), [S * S % P])
        for m in F.roots(cub):
            if (m - x3) % P == 0:
                continue
            w = S * ((x2 - m) % P) % P * pow((m - x3) % P, -1, P) % P
            out.append({12186: m, 24908: (y1 + w) % P})
    return out


base = report(v0, 'start (%s)' % src)[0]
best, bestv = base, list(v0)
t0 = time.time()
for pair in [('x1', 'y1'), ('x2', 'y2'), ('x3', 'y3'), ('x2', 'y1'), ('x1', 'y2')]:
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
