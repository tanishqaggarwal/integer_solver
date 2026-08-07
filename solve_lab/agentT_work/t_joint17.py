#!/usr/bin/env python3
"""AUDIT T32c -- the |S|=17 JOINT two-wire solve.

T32b named the block: on the only wire pair the residue admits, 398 of 400 sampled bivariate roots
clear ((x10261-x8912)-(13040669*x27539)) by direct recomputation and EVERY one of them breaks the
same two atoms, so the net count goes 3 -> 4.  So solve the three conditions TOGETHER: find
(t_w,t_v) zeroing the target's fitted polynomial mod c_target AND leaving both collateral atoms
divisible by their own c.  Exhaustive per prime power (iterate t_w over [0,q^e), root-find in t_v,
intersect), then CRT, then the usual two guards: direct recomputation and the global nonzero count.
"""
import os, sys, json, itertools, random, collections
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
import t_poly as TP
E = C.E; SL = C.SL; SHIFT = C.SHIFT; p = C.p; NV = C.NV
relift = C.relift; vars_of = C.vars_of; atomvalvars = C.atomvalvars
influences = C.influences; nzcount = C.nzcount; factor = C.factor; crt_list = C.crt_list

vv = [0]*NV
for k, val in json.load(open(os.path.join(T, 'close_T17w.json'))).items():
    vv[int(k[2:])] = int(val)
relift(vv)
BASE = nzcount(vv)
VV0 = vv[:]
r = E.run(vv)
ZERO0 = set(i for i, x in enumerate(r) if x == 0)
print('|S|=17 two-wire end state: global nonzero %d' % BASE)

TGTW = ('x24468', 'x18956')
A = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0
     and not any(t in a for t in TGTW)][0]
COLL = ['((x30973*x25513)-(1707229*x4242))', '((x38170*x15286)-(2264251*x9216))']
GROUP = [A] + COLL
W, V = 10261, 27156
print('wire pair (x%d, x%d)   [the ONLY pair the residue admits -- 2 wires]' % (W, V))
for a in GROUP:
    c = abs(SL[a])//p
    print('   %-52s c=%-10d = %s  current residual %s' % (
        a[:52], c, '*'.join('%d^%d' % qe for qe in sorted(factor(c).items())),
        'NONZERO' if r[E.residx[a]] else 'zero'))

# ------------------------------------------------------------------ fit all three
CF = {}
for a in GROUP:
    cf = C.fit2(vv, E.residx[a], W, V)
    if cf is None:
        print('   !! %s not p-divisible on the grid' % a[:40]); sys.exit(1)
    ok = all(C.probe2(vv, E.residx[a], W, tw, V, tv)//p == C.peval2_exact(cf, tw, tv)
             for tw, tv in ((7, 11), (1234, 99), (55555, 4), (3, 91237)))
    cross = any(cf[k][l] for k in range(1, len(cf)) for l in range(1, len(cf)))
    print('   fit %-40s validated-vs-recomputation %s  cross-term %s'
          % (a[:40], 'OK' if ok else 'MISMATCH', cross))
    if not ok:
        sys.exit(1)
    CF[a] = cf

# ------------------------------------------------------------------ joint root sets
PP = {}
for a in GROUP:
    for q, e in factor(abs(SL[a])//p).items():
        PP[q] = max(PP.get(q, 0), e)
print('\nprime powers to satisfy: %s' % sorted(PP.items()))

def binom_mod(t, k, m):
    b = 1
    for j in range(k):
        b = b*(t-j)//(j+1)
    return b % m

def rootsets(q, e):
    """EXHAUSTIVE over t_w in [0,q^e): all (t_w,t_v) satisfying every applicable condition."""
    m = q**e
    need = [a for a in GROUP if (abs(SL[a])//p) % q == 0]
    mods = {a: q**min(e, factor(abs(SL[a])//p)[q]) for a in need}
    D = len(CF[GROUP[0]])-1
    out = []
    for tw in range(m):
        Bw = [binom_mod(tw, k, m) for k in range(D+1)]
        cand = None
        for a in need:
            ma = mods[a]
            u = [sum(CF[a][k][l]*Bw[k] for k in range(D+1)) % ma for l in range(D+1)]
            mo = TP.newton_to_mono(u, ma) if ma > D else None
            if mo is None:
                rs = set(b for b in range(ma)
                         if sum(u[l]*binom_mod(b, l, ma) for l in range(D+1)) % ma == 0)
            else:
                rr = TP.roots_pp(mo, ma, q, min(e, factor(abs(SL[a])//p)[q]))
                rs = set(range(ma)) if rr == 'ALL' else set(rr)
            rs = set(b for b in range(m) if b % ma in rs) if ma < m else rs
            cand = rs if cand is None else (cand & rs)
            if not cand:
                break
        if cand:
            out.extend((tw, b) for b in cand)
    return out

per = []
for q in sorted(PP):
    e = PP[q]
    rs = rootsets(q, e)
    print('  mod %-8d : %d joint (t_w,t_v) pairs  [exhaustive over t_w]' % (q**e, len(rs)))
    per.append((q**e, rs))
    if not rs:
        print('\n>>> NO JOINT ROOT mod %d.  The three conditions are simultaneously unsatisfiable'
              ' on this wire pair, exhaustively.' % q**e)
        sys.exit(0)

n = 1
for _, rs in per:
    n *= len(rs)
print('\njoint solution count mod lcm = %d' % n)

# ------------------------------------------------------------------ guards
rnd = random.Random(77)
mods = [m for m, _ in per]; sets = [rs for _, rs in per]
seen = set(); tried = cleared = 0; best = None; pat = collections.Counter(); WIN = None
while tried < 600:
    pick = [rs[rnd.randrange(len(rs))] for rs in sets]
    tw = crt_list([(x, m) for (x, _), m in zip(pick, mods)])
    tv = crt_list([(y, m) for (_, y), m in zip(pick, mods)])
    if tw is None or tv is None or (tw, tv) in seen:
        tried += 1; continue
    seen.add((tw, tv)); tried += 1
    if (tw, tv) == (0, 0):
        continue
    good = all(C.probe2(vv, E.residx[a], W, tw, V, tv) % abs(SL[a]) == 0 for a in GROUP)
    if not good:
        continue
    cleared += 1
    vv[W] += p*tw; vv[V] += p*tv
    nn = nzcount(vv)
    rn = E.run(vv)
    new = tuple(sorted(i for i in ZERO0 if rn[i]))
    if nn < BASE:
        print('\n*** x%d += p*%d , x%d += p*%d  -> global %d -> %d' % (W, tw, V, tv, BASE, nn))
        WIN = (tw, tv, nn)
        break
    vv[:] = VV0
    pat[(new, nn)] += 1
    if best is None or nn < best:
        best = nn
print('\ntried %d CRT combinations, %d verified by direct recomputation on all three atoms'
      % (tried, cleared))
if WIN:
    relift(vv)
    rr = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(rr) if x]
    print('|S|=17 now has %d nonzero atoms:' % len(nz))
    for a in nz:
        print('   ', a[:110])
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
              open(os.path.join(T, 'close_T17j.json'), 'w'))
    print('dumped close_T17j.json -- run checker.py on it')
else:
    print('best global %s (base %d); collateral patterns:' % (best, BASE))
    for (new, nn) in list(pat)[:5]:
        print('   ->%d  newly nonzero %d: %s' % (nn, len(new), [E.res[i][:64] for i in new[:4]]))
