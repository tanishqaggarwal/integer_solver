#!/usr/bin/env python3
"""AUDIT T32 -- closeS4 + a TWO-WIRE fallback.

closeS4 stops when no SINGLE-wire shift strictly decreases the global nonzero-atom count.  T31
showed that at |S|=8 the surviving condition has no univariate root on any of its 6 wires but a
bivariate root on 10 of its 15 wire pairs, and that one of them clears it with zero collateral
(checker-verified 39,018).  So: run the same closure loop with a two-wire pass that fires whenever
the single-wire pass stalls, and take it to |S| = 17, the campaign's decisive open case.

Guards, unchanged from closeS4 / L / P:
  * every candidate root is verified by DIRECT RECOMPUTATION of the atom (probe, not the fit);
  * a shift is accepted only if the TOTAL nonzero-atom count strictly decreases;
  * the fitted bivariate polynomial is validated against direct recomputation before it is used.
Usage:  python3 t_close2w.py <tag> <n|explicit,comma,list>
"""
import os, sys, json, time, itertools, random, collections
from math import gcd
T = '/home/user/integer_solver/solve_lab/agentT_work'
L = os.path.join(T, 'mirror', 'L')
sys.path.insert(0, T)
import t_poly as TP
os.chdir(L); sys.path.insert(0, L)
g = {'__name__': 'drv'}
exec(compile(open(os.path.join(L, 'closeS4.py')).read().split("if __name__")[0], 'c4', 'exec'), g)
E = g['E']; SL = g['SL']; SHIFT = g['SHIFT']; p = g['p']; NV = g['NV']; M = g['M']
relift = g['relift']; vars_of = g['vars_of']; atomvalvars = g['atomvalvars']
influences = g['influences']; nzcount = g['nzcount']; assignment = g['assignment']
ORIENT = g['ORIENT']; T1 = g['T1']; T2 = g['T2']
factor = g['factor']; crt_list = g['crt_list']; fitc = g['fitc']; roots_c = g['roots_c']
solve_group3 = g['solve_group3']

DEG = 5
SMALLM = 400          # prime powers up to this are enumerated exhaustively as q^2 pairs
SAMP = 120            # t_w samples per prime power above SMALLM
GUARD = 250           # global-guard evaluations per (atom,pair) at most
rnd = random.Random(20260807)

# ------------------------------------------------------------------ bivariate machinery
def probe2(vv, i, w, tw, v, tv):
    ow, ov = vv[w], vv[v]
    vv[w] = ow + p*tw; vv[v] = ov + p*tv
    y = E.run(vv)[i]
    vv[w] = ow; vv[v] = ov
    return y

def fit2(vv, i, w, v, D=DEG):
    ys = []
    for a_ in range(D+1):
        row = []
        for b_ in range(D+1):
            y = probe2(vv, i, w, a_, v, b_)
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

def binom_mod(t, k, m):
    b = 1
    for j in range(k):
        b = b*(t-j)//(j+1)
    return b % m

def collapse(cf, tw, m):
    """Newton coeffs of the univariate polynomial in t_v obtained by fixing t_w."""
    D = len(cf)-1
    B = [binom_mod(tw, k, m) for k in range(D+1)]
    return [sum(cf[k][l]*B[k] for k in range(D+1)) % m for l in range(D+1)]

def pairroots(cf, q, e, want):
    """(t_w,t_v) pairs that zero the fitted poly mod q^e.  Exhaustive for small q^e, sampled
    above it.  Second return value says whether the search was exhaustive."""
    m = q**e
    D = len(cf)-1
    if m <= SMALLM:
        B = [[binom_mod(t, k, m) for k in range(D+1)] for t in range(m)]
        out = []
        for a_ in range(m):
            u = [sum(cf[k][l]*B[a_][k] for k in range(D+1)) % m for l in range(D+1)]
            if not any(u):
                out.extend((a_, b_) for b_ in range(m)); continue
            for b_ in range(m):
                if sum(u[l]*B[b_][l] for l in range(D+1)) % m == 0:
                    out.append((a_, b_))
        return out, True
    out = []; tried = 0
    while len(out) < want and tried < SAMP:
        tried += 1
        a_ = rnd.randrange(m)
        u = collapse(cf, a_, m)
        mo = TP.newton_to_mono(u, m)
        if mo is None:
            continue
        rs = TP.roots_pp(mo, m, q, e)
        if rs == 'ALL':
            out.extend((a_, rnd.randrange(m)) for _ in range(4))
        else:
            out.extend((a_, b_) for b_ in rs)
    return out, False

# ------------------------------------------------------------------ the two-wire pass
def two_wire(vv, viol, base, log):
    """try to clear ONE violated atom with a two-wire shift; returns True on an accepted shift."""
    for a in viol:
        i = E.residx[a]; c = abs(SL[a])//p
        if c <= 1:
            continue
        WS = sorted(set(x for x in vars_of(E.atoms[a]) if x in SHIFT) |
                    set(x for x in atomvalvars[a] if x in SHIFT))
        WS = [w for w in WS if influences(vv, a, w)]
        FC = sorted(factor(c).items())
        log('  atom %s  c=%d=%s  wires=%d pairs=%d' % (
            a[:56], c, '*'.join('%d^%d' % qe for qe in FC), len(WS), len(WS)*(len(WS)-1)//2))
        for w, v in itertools.combinations(WS, 2):
            cf = fit2(vv, i, w, v)
            if cf is None:
                continue
            ok = all(probe2(vv, i, w, tw, v, tv)//p == peval2_exact(cf, tw, tv)
                     for tw, tv in ((7, 11), (1234, 99)))
            if not ok:
                log('     (x%d,x%d) FIT DISAGREES WITH RECOMPUTATION -- skipped' % (w, v))
                continue
            per = []; dead = None
            for q, e in FC:
                rs, ex = pairroots(cf, q, e, 24)
                per.append((q**e, rs, ex))
                if not rs:
                    dead = (q, ex); break
            if dead:
                log('     (x%-6d,x%-6d) no root mod %d %s' % (
                    w, v, dead[0], '(exhaustive)' if dead[1] else '(sampled)'))
                continue
            mods = [m for m, _, _ in per]; sets = [rs for _, rs, _ in per]
            seen = set(); tried = cleared = 0
            while tried < GUARD:
                pick = [rs[rnd.randrange(len(rs))] for rs in sets]
                tw = crt_list([(x, m) for (x, _), m in zip(pick, mods)])
                tv = crt_list([(y, m) for (_, y), m in zip(pick, mods)])
                if tw is None or tv is None or (tw, tv) in seen:
                    tried += 1; continue
                seen.add((tw, tv)); tried += 1
                if (tw, tv) == (0, 0):
                    continue
                if probe2(vv, i, w, tw, v, tv) % (c*p) != 0:
                    continue
                cleared += 1
                snap = vv[:]
                vv[w] += p*tw; vv[v] += p*tv
                n = nzcount(vv)
                if n < base:
                    log('     *** (x%d += p*%d , x%d += p*%d) CLEARS, global %d -> %d'
                        % (w, tw, v, tv, base, n))
                    return True
                vv[:] = snap
            log('     (x%-6d,x%-6d) %d roots tried, %d verified, none passed the global guard'
                % (w, v, tried, cleared))
    return False

# ------------------------------------------------------------------ the closure loop
def close2w(S, tag, outer_max=14, logf=None):
    def log(s):
        print(s, flush=True)
        if logf:
            logf.write(s+'\n'); logf.flush()
    v, isl, valn = assignment(set(S), ORIENT); v[24468] = T1; v[18956] = T2
    vv = [0]*NV
    for k, x in v.items():
        vv[k] = x
    for rd in range(60):                       # closeS4's greedy fixpoint, verbatim
        bad = relift(vv)
        if not bad:
            break
        r = E.run(vv); fx = 0
        for a in bad:
            i = E.residx[a]; cur = r[i]; sm = abs(SL[a])
            if cur % p:
                continue
            imm = [q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old = vv[w]; vv[w] = old+p; d = E.run(vv)[i]-cur; vv[w] = old
                if d == 0:
                    continue
                gg = gcd(d, sm)
                if cur % gg:
                    continue
                mm = sm//gg
                t = (-(cur//gg))*pow((d//gg) % mm, -1, mm) % mm if mm > 1 else 0
                vv[w] = old+p*t; fx += 1; break
        if fx == 0:
            break
    gen = 0
    for outer in range(outer_max):
        base = nzcount(vv); r = E.run(vv); gen += 1
        viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0]
        log('outer %d: global nonzero %d, violated c-conditions %d' % (outer, base, len(viol)))
        if not viol:
            break
        wires = collections.defaultdict(list)
        for a in viol:
            for w in (set(q for q in vars_of(E.atoms[a]) if q in SHIFT) |
                      set(q for q in atomvalvars[a] if q in SHIFT)):
                wires[w].append(a)
        prog = 0
        for w, ats in sorted(wires.items(), key=lambda kv: -len(kv[1])):
            V = [a for a in ats if influences(vv, a, w)]
            if not V:
                continue
            t = solve_group3(vv, V, w, gen, base)
            if t:
                prog += 1; base = nzcount(vv); gen += 1
        if prog:
            log('   single-wire pass: %d accepted shift(s)' % prog)
            continue
        log('   single-wire pass STALLED -> two-wire pass')
        r = E.run(vv)
        viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0]
        if not two_wire(vv, viol, base, log):
            log('   two-wire pass also stalled -> stop')
            break
        gen += 1
    relift(vv); r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
              open(os.path.join(T, 'close_%s.json' % tag), 'w'))
    return nz

if __name__ == '__main__':
    tag = sys.argv[1]; spec = sys.argv[2]
    if ',' in spec:
        S = [int(x) for x in spec.split(',')]
    else:
        n = int(spec)
        r7 = random.Random(7)
        S = [24601, 2081] if n == 2 else r7.sample(M['live'], n)
    lf = open(os.path.join(T, 't_close2w_%s.log' % tag), 'w')
    lf.write('S = %s\n' % S)
    print('S = %s' % S, flush=True)
    t0 = time.time()
    nz = close2w(S, tag, logf=lf)
    el = time.time()-t0
    msg = ('|S|=%-3d %-6s  NONZERO ATOMS = %d of %d   WALL = %.1f s  -> close_%s.json'
           % (len(S), tag, len(nz), len(E.res), el, tag))
    print(msg, flush=True); lf.write(msg+'\n')
    for a in nz[:10]:
        print('    ', a[:110], flush=True); lf.write('     '+a[:110]+'\n')
    lf.close()
