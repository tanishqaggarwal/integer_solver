#!/usr/bin/env python3
"""AUDIT T33 -- the general closer: closeS4 + a JOINT two-wire pass.

T31: at |S|=8 the residue had no univariate root but a bivariate one with zero collateral.
T32: at |S|=17 the residue had bivariate roots on its only pair, all blocked by exactly two
     collateral atoms; solving those three conditions TOGETHER on the same pair closed it
     (checker-verified 39,018).
This makes that a method rather than a one-off: whenever a two-wire root is rejected by the global
guard, the atoms it broke are folded into the system and the pair is re-solved.  Guards unchanged:
direct recomputation of every atom in the group, then the global nonzero-atom count must drop.

Usage:  python3 t_close2wj.py <tag> <n | comma,list>
"""
import os, sys, json, time, itertools, random, collections
from math import gcd
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
import t_poly as TP
E = C.E; SL = C.SL; SHIFT = C.SHIFT; p = C.p; NV = C.NV; M = C.M
relift = C.relift; vars_of = C.vars_of; atomvalvars = C.atomvalvars
influences = C.influences; nzcount = C.nzcount; factor = C.factor; crt_list = C.crt_list
assignment = C.assignment; ORIENT = C.ORIENT; T1 = C.T1; T2 = C.T2
solve_group3 = C.solve_group3

EXCAP = 60000        # prime powers up to this are solved EXHAUSTIVELY over t_w
SAMPW = 400          # t_w samples above it
TRIES = 200          # CRT combinations tested per group
ROUNDS = 4           # collateral-refinement rounds per pair
rnd = random.Random(20260807)

def binom_mod(t, k, m):
    b = 1
    for j in range(k):
        b = b*(t-j)//(j+1)
    return b % m

def tv_roots(cf, tw, ma, q, ee, D):
    u = [sum(cf[k][l]*binom_mod(tw, k, ma) for k in range(D+1)) % ma for l in range(D+1)]
    if ma > D:
        mo = TP.newton_to_mono(u, ma)
        if mo is not None:
            rr = TP.roots_pp(mo, ma, q, ee)
            return set(range(ma)) if rr == 'ALL' else set(rr)
    return set(b for b in range(ma)
               if sum(u[l]*binom_mod(b, l, ma) for l in range(D+1)) % ma == 0)

def joint_rootsets(CF, GROUP, q, e):
    """(t_w,t_v) pairs mod q^e satisfying every applicable condition; exhaustive when it can be."""
    m = q**e
    need = [a for a in GROUP if (abs(SL[a])//p) % q == 0]
    if not need:
        return [(0, 0)], True, m           # unconstrained at this prime: t=0 is safe
    ex = {a: factor(abs(SL[a])//p)[q] for a in need}
    D = len(CF[GROUP[0]])-1
    exhaustive = m <= EXCAP
    tws = range(m) if exhaustive else [rnd.randrange(m) for _ in range(SAMPW)]
    out = []
    for tw in tws:
        cand = None
        for a in need:
            ma = q**min(e, ex[a])
            rs = tv_roots(CF[a], tw, ma, q, min(e, ex[a]), D)
            if ma < m:
                rs = set(b for b in range(m) if b % ma in rs)
            cand = rs if cand is None else (cand & rs)
            if not cand:
                break
        if cand:
            out.extend((tw, b) for b in cand)
            if not exhaustive and len(out) > 40:
                break
    return out, exhaustive, m

def joint_pair(vv, a0, base, log):
    """try to clear atom a0 with a joint two-wire shift; True if an accepted shift was applied."""
    WS = sorted(set(x for x in vars_of(E.atoms[a0]) if x in SHIFT) |
                set(x for x in atomvalvars[a0] if x in SHIFT))
    WS = [w for w in WS if influences(vv, a0, w)]
    c0 = abs(SL[a0])//p
    log('  atom %s  c=%d  wires=%d pairs=%d' % (a0[:56], c0, len(WS), len(WS)*(len(WS)-1)//2))
    r0 = E.run(vv)
    ZERO0 = set(i for i, x in enumerate(r0) if x == 0)
    for w, v in itertools.combinations(WS, 2):
        GROUP = [a0]; CF = {}
        for rd in range(ROUNDS):
            okfit = True
            for a in GROUP:
                if a in CF:
                    continue
                cf = C.fit2(vv, E.residx[a], w, v)
                if cf is None or not all(
                        C.probe2(vv, E.residx[a], w, tw, v, tv)//p == C.peval2_exact(cf, tw, tv)
                        for tw, tv in ((7, 11), (1234, 99), (55555, 4))):
                    okfit = False; break
                CF[a] = cf
            if not okfit:
                log('     (x%d,x%d) no valid fit -- skipped' % (w, v)); break
            PP = {}
            for a in GROUP:
                for q, e in factor(abs(SL[a])//p).items():
                    PP[q] = max(PP.get(q, 0), e)
            per = []; dead = None
            for q in sorted(PP):
                rs, exh, m = joint_rootsets(CF, GROUP, q, PP[q])
                per.append((m, rs))
                if not rs:
                    dead = (m, exh); break
            if dead:
                log('     (x%-6d,x%-6d) |group|=%d  NO JOINT ROOT mod %d %s'
                    % (w, v, len(GROUP), dead[0], '(exhaustive)' if dead[1] else '(sampled)'))
                break
            mods = [m for m, _ in per]; sets = [rs for _, rs in per]
            seen = set(); tried = cleared = 0; newatoms = collections.Counter()
            while tried < TRIES:
                pick = [rs[rnd.randrange(len(rs))] for rs in sets]
                tw = crt_list([(x, m) for (x, _), m in zip(pick, mods)])
                tv = crt_list([(y, m) for (_, y), m in zip(pick, mods)])
                if tw is None or tv is None or (tw, tv) in seen:
                    tried += 1; continue
                seen.add((tw, tv)); tried += 1
                if (tw, tv) == (0, 0):
                    continue
                if not all(C.probe2(vv, E.residx[a], w, tw, v, tv) % abs(SL[a]) == 0
                           for a in GROUP):
                    continue
                cleared += 1
                snap = vv[:]
                vv[w] += p*tw; vv[v] += p*tv
                n = nzcount(vv)
                if n < base:
                    log('     *** (x%d += p*%d , x%d += p*%d) |group|=%d  global %d -> %d'
                        % (w, tw, v, tv, len(GROUP), base, n))
                    return True
                rn = E.run(vv)
                for i in ZERO0:
                    if rn[i]:
                        newatoms[E.res[i]] += 1
                vv[:] = snap
            add = [a for a, _ in newatoms.most_common() if a not in GROUP and a in SL and SL[a]
                   and abs(SL[a])//p > 1][:3]
            log('     (x%-6d,x%-6d) |group|=%d  %d roots verified, guard rejected all; '
                'collateral %s' % (w, v, len(GROUP), cleared,
                                   [a[:44] for a in add] or 'none addable'))
            if not add:
                break
            GROUP += add
    return False

def handleless_pass(vv, base, log):
    """T33c: atoms with NO handle cannot absorb anything and must be EXACTLY zero, so closeS4's
    condition machinery (which is indexed by SL, i.e. by handle) never touches them.  At |S|=32
    two such atoms are nonzero from construction and each has one admitted wire on which R(t) is
    LINEAR with an exact integer root.  Zero them, with the usual global guard."""
    r = E.run(vv)
    HL = [a for a in E.res if r[E.residx[a]] and a not in SL]
    if not HL:
        return 0
    log('   handle-less pass: %d nonzero handle-less atom(s)' % len(HL))
    fixed = 0
    for a in HL:
        i = E.residx[a]
        ws = sorted(set(x for x in vars_of(E.atoms[a]) if x in SHIFT) |
                    set(x for x in atomvalvars[a] if x in SHIFT))
        ws = [w for w in ws if influences(vv, a, w)]
        for w in ws:
            ys = []
            for t in range(6):
                ow = vv[w]; vv[w] = ow + p*t; ys.append(E.run(vv)[i]); vv[w] = ow
            d = [ys[:]]
            for k in range(5):
                d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
            cf = [d[k][0] for k in range(6)]
            top = max([k for k in range(6) if cf[k]], default=0)
            cand = []
            if top == 1 and cf[1] and cf[0] % cf[1] == 0:
                cand = [-cf[0]//cf[1]]
            elif top >= 2:
                log('      x%d: R(t) degree %d -- integer root not attempted' % (w, top))
            for t in cand:
                snap = vv[:]
                vv[w] = vv[w] + p*t
                if E.run(vv)[i] == 0 and nzcount(vv) < base:
                    log('      *** x%d += p*%d zeroes %s (exact, recomputed)' % (w, t, a[:44]))
                    fixed += 1; base = nzcount(vv); break
                vv[:] = snap
            if fixed:
                break
    return fixed


def close(S, tag, outer_max=16, logf=None):
    def log(s):
        print(s, flush=True)
        if logf:
            logf.write(s+'\n'); logf.flush()
    v, isl, valn = assignment(set(S), ORIENT); v[24468] = T1; v[18956] = T2
    vv = [0]*NV
    for k, x in v.items():
        vv[k] = x
    for rd in range(60):
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
    TGT = ('x24468', 'x18956')
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
            log('   single-wire pass: %d accepted' % prog); continue
        if handleless_pass(vv, nzcount(vv), log):
            continue
        log('   single-wire STALLED -> joint two-wire pass')
        r = E.run(vv)
        viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0
                and not any(t in a for t in TGT)]
        if not viol:
            log('   only the two TARGET congruences remain -- CLOSED'); break
        if not any(joint_pair(vv, a, base, log) for a in viol):
            log('   joint two-wire pass also stalled -> stop'); break
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
        n = int(spec); r7 = random.Random(7)
        S = [24601, 2081] if n == 2 else r7.sample(M['live'], n)
    lf = open(os.path.join(T, 't_close2wj_%s.log' % tag), 'w')
    lf.write('S = %s\n' % S); print('S = %s' % S, flush=True)
    t0 = time.time(); nz = close(S, tag, logf=lf); el = time.time()-t0
    msg = ('|S|=%-3d %-6s  NONZERO ATOMS = %d of %d   WALL = %.1f s  -> close_%s.json'
           % (len(S), tag, len(nz), len(E.res), el, tag))
    print(msg, flush=True); lf.write(msg+'\n')
    for a in nz[:10]:
        print('    ', a[:110], flush=True); lf.write('     '+a[:110]+'\n')
    lf.close()
