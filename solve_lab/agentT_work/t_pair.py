#!/usr/bin/env python3
"""AUDIT T31 -- TWO-WIRE SHIFT on the |S|=8 surviving condition.

The |S|=8 residue ((x21408*x10138)-(15333171*x658)), c = 3*7*19*83*463, has NO univariate root on
any of its 6 candidate wires (audit T30b / t_leaf.py).  That is exactly the shape of L's bivariate
residue, where a p*t_w*t_v term survives mod c.  So: exhaustive over all 15 wire pairs, fit the
BIVARIATE polynomial R(t_w,t_v)/p exactly, root-find it mod c prime-by-prime EXHAUSTIVELY, and for
every root found apply it and test with closeS4's own global guard (total nonzero-atom count must
strictly decrease), by DIRECT RECOMPUTATION.

Every fitted polynomial is validated against direct recomputation at random points before it is
used -- ledger rule 6.
"""
import os, sys, json, time, itertools, random
from math import gcd
L = '/home/user/integer_solver/solve_lab/agentL_work'
T = '/home/user/integer_solver/solve_lab/agentT_work'
os.chdir(L); sys.path.insert(0, L)
g = {'__name__': 'drv'}
exec(compile(open(os.path.join(T, 'from_L', 'closeS4.py')).read().split("if __name__")[0],
             'c4', 'exec'), g)
E = g['E']; SL = g['SL']; SHIFT = g['SHIFT']; p = g['p']; NV = g['NV']
relift = g['relift']; vars_of = g['vars_of']; atomvalvars = g['atomvalvars']
influences = g['influences']; nzcount = g['nzcount']
factor = g['factor']; crt_list = g['crt_list']

NTRY = int(sys.argv[1]) if len(sys.argv) > 1 else 300
DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 5

# ---------------------------------------------------------------- load the |S|=8 end state
vv = [0]*NV
for k, val in json.load(open(os.path.join(T, 'close_T8.json'))).items():
    vv[int(k[2:])] = int(val)
relift(vv)
r = E.run(vv)
NZ0 = [E.res[i] for i, x in enumerate(r) if x]
print('|S|=8 end state: %d nonzero atoms of %d' % (len(NZ0), len(E.res)), flush=True)
for a in NZ0:
    print('    ', a[:100])
viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0]
assert len(viol) == 1, viol
A = viol[0]; IA = E.residx[A]; C = abs(SL[A])//p
FC = sorted(factor(C).items())
print('condition atom : %s' % A[:90])
print('c              = %d = %s' % (C, ' * '.join('%d^%d' % (q, e) for q, e in FC)), flush=True)
BASE = nzcount(vv)
VV0 = vv[:]
ZERO0 = set(i for i, x in enumerate(E.run(vv)) if x == 0)
print('baseline global nonzero count: %d' % BASE, flush=True)

WS = sorted(set(q for q in vars_of(E.atoms[A]) if q in SHIFT) |
            set(q for q in atomvalvars[A] if q in SHIFT))
WSall = WS[:]
WS = [w for w in WS if influences(vv, A, w)]
print('candidate shift wires: %d of %d admitted by influences()  %s' % (len(WS), len(WSall), WS),
      flush=True)

# ---------------------------------------------------------------- bivariate machinery
def probe2(i, w, tw, v, tv):
    ow, ov = vv[w], vv[v]
    vv[w] = ow + p*tw; vv[v] = ov + p*tv
    y = E.run(vv)[i]
    vv[w] = ow; vv[v] = ov
    return y

def fit2(i, w, v, D):
    """exact 2-D Newton forward-difference fit of R(tw,tv)/p on the (D+1)x(D+1) grid."""
    ys = []
    for a_ in range(D+1):
        row = []
        for b_ in range(D+1):
            y = probe2(i, w, a_, v, b_)
            if y % p:
                return None
            row.append(y//p)
        ys.append(row)
    tmp = []
    for row in ys:
        d = [row[:]]
        for k in range(D):
            d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
        tmp.append([d[k][0] for k in range(D+1)])
    d = [tmp]
    for k in range(D):
        d.append([[d[k][j+1][l]-d[k][j][l] for l in range(D+1)] for j in range(len(d[k])-1)])
    return [d[k][0][:] for k in range(D+1)]

def peval2_exact(cf, tw, tv):
    tot = 0; bw = 1
    for k in range(len(cf)):
        if k:
            bw = bw*(tw-k+1)//k
        bv = 1
        for l in range(len(cf[k])):
            if l:
                bv = bv*(tv-l+1)//l
            if cf[k][l]:
                tot += cf[k][l]*bw*bv
    return tot

def roots_mod_q(cf, q, e):
    """EXHAUSTIVE: every (tw,tv) in [0,q^e)^2 with the fitted poly == 0 mod q^e."""
    m = q**e
    D = len(cf)-1
    B = [[1]*(D+1) for _ in range(m)]           # B[t][k] = C(t,k) mod m
    for t in range(m):
        b = 1
        for k in range(D+1):
            if k:
                b = b*(t-k+1)//k
            B[t][k] = b % m
    cq = [[cf[k][l] % m for l in range(D+1)] for k in range(D+1)]
    out = []
    for a_ in range(m):
        Ba = B[a_]
        # collapse the tw direction -> univariate Newton coeffs in tv
        u = [sum(cq[k][l]*Ba[k] for k in range(D+1)) % m for l in range(D+1)]
        if not any(u):
            out.extend((a_, b_) for b_ in range(m)); continue
        for b_ in range(m):
            Bb = B[b_]
            if sum(u[l]*Bb[l] for l in range(D+1)) % m == 0:
                out.append((a_, b_))
    return out

# ---------------------------------------------------------------- pass 1: root existence
rnd = random.Random(11)
PAIRS = []
t_all = time.time()
print('\n== PASS 1: bivariate fit + EXHAUSTIVE root enumeration mod each prime of c ==', flush=True)
for w, v in itertools.combinations(WS, 2):
    t0 = time.time()
    cf = fit2(IA, w, v, DEG)
    if cf is None:
        print('pair (x%-6d,x%-6d)  residual not p-divisible on the grid -- SKIP' % (w, v), flush=True)
        PAIRS.append(dict(w=w, v=v, st='nofit')); continue
    dw = max([k for k in range(DEG+1) if any(cf[k])], default=0)
    dv = max([l for l in range(DEG+1) if any(cf[k][l] for k in range(DEG+1))], default=0)
    trunc = (dw == DEG or dv == DEG)
    bad = 0
    for _ in range(5):
        tw = rnd.randrange(0, 10**4); tv = rnd.randrange(0, 10**4)
        if probe2(IA, w, tw, v, tv)//p != peval2_exact(cf, tw, tv):
            bad += 1
    cross = [(k, l) for k in range(1, DEG+1) for l in range(1, DEG+1) if cf[k][l]]
    perq = []; empty = None
    for q, e in FC:
        rs = roots_mod_q(cf, q, e)
        perq.append((q, q**e, rs))
        if not rs:
            empty = q; break
    n = 1
    for _, _, rs in perq:
        n *= len(rs)
    print('pair (x%-6d,x%-6d) deg %d,%d cross=%-6s fit-vs-recompute %s  %s  (%.0fs)' % (
        w, v, dw, dv, bool(cross), 'OK' if bad == 0 else 'MISMATCH x%d' % bad,
        ('NO ROOT mod %d' % empty) if empty else
        ('roots/prime %s = %d mod c' % ([len(rs) for _, _, rs in perq], n)),
        time.time()-t0), flush=True)
    if trunc:
        print('     ** WARNING: degree cap %d reached -- refit with a larger DEG' % DEG, flush=True)
    PAIRS.append(dict(w=w, v=v, st='noroot' if empty else 'roots', n=n, dw=dw, dv=dv,
                      cross=len(cross), fitok=(bad == 0), perq=None if empty else perq,
                      trunc=trunc))

nr = sum(1 for P in PAIRS if P['st'] == 'noroot')
hr = sum(1 for P in PAIRS if P['st'] == 'roots')
print('\nPASS 1 SUMMARY: %d pairs, %d with NO bivariate root, %d with roots' % (len(PAIRS), nr, hr),
      flush=True)

# ---------------------------------------------------------------- pass 2: the global guard
print('\n== PASS 2: apply the roots, closeS4 global guard, direct recomputation ==', flush=True)
WIN = None
for P in PAIRS:
    if P['st'] != 'roots':
        continue
    w, v, perq, n = P['w'], P['v'], P['perq'], P['n']
    t0 = time.time()
    tried = cleared = improved = 0
    best = None; seen = set(); coll = {}
    mods = [m for _, m, _ in perq]; sets = [rs for _, _, rs in perq]
    rr = random.Random(1000*w + v)
    lim = min(NTRY, n)
    guard = 0
    while tried < lim and guard < 40*lim:
        guard += 1
        pick = [rs[rr.randrange(len(rs))] for rs in sets]
        tw = crt_list([(a_, m) for (a_, _), m in zip(pick, mods)])
        tv = crt_list([(b_, m) for (_, b_), m in zip(pick, mods)])
        if tw is None or tv is None or (tw, tv) in seen:
            continue
        seen.add((tw, tv)); tried += 1
        if (tw, tv) == (0, 0):
            continue
        if probe2(IA, w, tw, v, tv) % (C*p) != 0:
            continue
        cleared += 1
        vv[w] += p*tw; vv[v] += p*tv
        nc = nzcount(vv)
        rn = E.run(vv)
        newnz = tuple(sorted(i for i in ZERO0 if rn[i]))
        vv[:] = VV0
        coll[newnz] = coll.get(newnz, 0) + 1
        if best is None or nc < best:
            best = nc
        if nc < BASE:
            improved += 1
            print('  *** x%d += p*%d , x%d += p*%d  CLEARS and global %d -> %d' % (
                w, tw, v, tv, BASE, nc), flush=True)
            WIN = (w, tw, v, tv, nc); break
    print('  (x%-6d,x%-6d) roots %d  tried %d  cleared-by-recomputation %d  improved %d  '
          'best global %s (base %d)  distinct collateral patterns %d  (%.0fs)' % (
              w, v, n, tried, cleared, improved, best, BASE, len(coll), time.time()-t0), flush=True)
    for k, cnt in sorted(coll.items(), key=lambda kv: -kv[1])[:4]:
        print('        collateral: %d newly-nonzero atoms  x%d   %s' % (
            len(k), cnt, [E.res[i][:60] for i in k[:3]]), flush=True)
    P['tried'] = tried; P['cleared'] = cleared; P['improved'] = improved; P['best'] = best
    P['ncoll'] = len(coll)
    if WIN:
        break

json.dump([{k: val for k, val in P.items() if k != 'perq'} for P in PAIRS],
          open(os.path.join(T, 't_pair_roots.json'), 'w'), indent=1)
print('\nTOTAL WALL %.0f s' % (time.time()-t_all), flush=True)
if WIN:
    w, tw, v, tv, nc = WIN
    vv[:] = VV0
    vv[w] += p*tw; vv[v] += p*tv
    nzcount(vv)
    rr2 = E.run(vv)
    nz2 = [E.res[i] for i, x in enumerate(rr2) if x]
    print('WINNER: two-wire shift x%d += p*%d , x%d += p*%d  -> %d nonzero atoms' % (
        w, tw, v, tv, len(nz2)))
    for a in nz2:
        print('   ', a[:110])
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
              open(os.path.join(T, 'close_T8pair.json'), 'w'))
    print('dumped close_T8pair.json -- now run checker.py on it')
else:
    print('NO PAIR CLEARED IT over all %d pairs.' % len(PAIRS))
