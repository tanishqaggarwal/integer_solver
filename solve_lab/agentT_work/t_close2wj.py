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
ROUNDS = 4           # collateral-refinement rounds per atom
PAIRCAP = 40         # wire pairs examined per round
OUTCAP  = 4000       # T39: max joint roots carried per prime power -- see joint_rootsets
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
    """(t_w,t_v) pairs mod q^e satisfying every applicable condition; exhaustive when it can be.

    T36 (rule 9 applied to my own solver): above EXCAP this SAMPLES t_w and root-finds t_v.  That
    is asymmetric, and it fails completely on the degenerate shape where the residue does not
    depend on t_v mod q -- then only ~1 in q values of t_w works and 400 random draws find it with
    probability 400/q.  That is what produced `NO JOINT ROOT mod 116507 (sampled)` at |S|=128.
    Fix: when the first orientation comes back empty, TRANSPOSE and scan the other way, which turns
    exactly that degenerate case into a direct univariate root-find.  A sampled empty result is now
    a statement about BOTH orientations."""
    m = q**e
    need = [a for a in GROUP if (abs(SL[a])//p) % q == 0]
    if not need:
        return [(0, 0)], True, m           # unconstrained at this prime: t=0 is safe
    ex = {a: factor(abs(SL[a])//p)[q] for a in need}
    D = len(CF[GROUP[0]])-1
    exhaustive = m <= EXCAP

    def scan(CFo, flip):
        tws = range(m) if exhaustive else [rnd.randrange(m) for _ in range(SAMPW)]
        out = []
        for tw in tws:
            cand = None
            for a in need:
                ma = q**min(e, ex[a])
                rs = tv_roots(CFo[a], tw, ma, q, min(e, ex[a]), D)
                if ma < m:
                    rs = set(b for b in range(m) if b % ma in rs)
                cand = rs if cand is None else (cand & rs)
                if not cand:
                    break
            if cand:
                # T39: OUTCAP.  When the collapsed polynomial vanishes IDENTICALLY mod q^e,
                # tv_roots returns all q^e residues; over an exhaustive t_w loop that is m^2
                # tuples -- up to 3.6e9 at EXCAP=60000.  This is what OOM-killed the |S|=192 run
                # at 14.6 GB RSS (and, earlier, the first |S|=128 run, losing its dump).  The
                # root SET is unchanged; only how many representatives we carry forward is
                # capped, and TRIES=200 of them are tested anyway.
                if len(cand) > OUTCAP:
                    cand = sorted(cand)[:OUTCAP]
                out.extend(((b, tw) if flip else (tw, b)) for b in cand)
                if len(out) > OUTCAP or (not exhaustive and len(out) > 40):
                    break
        return out

    out = scan(CF, False)
    if not out and not exhaustive:
        out = scan({a: transpose_cf(CF[a]) for a in need}, True)
    return out, exhaustive, m

# ------------------------------------------------------------------ T34: MIXED constraints
# An atom with NO handle cannot absorb anything, so its condition is R == 0 EXACTLY over Z, not
# R == 0 mod c.  A group can now carry both kinds on the same wire pair.  Observed shape at
# |S|=32: the exact condition depends on ONE of the two wires and is LINEAR there, so it PINS
# that wire to a single integer and the divisibility conditions are then univariate on the other.

def is_exact(a):
    """True if atom a has no handle -- it must be exactly zero over Z."""
    s = SL.get(a)
    return not (s and s % p == 0 and abs(s)//p > 1)

def newton_int_roots(u):
    """integer roots of the integer polynomial given by Newton coeffs u (sum u[k]*C(t,k)).
    'ALL' if identically zero, [] if none, None if the degree is beyond what is implemented."""
    u = list(u)
    top = max([k for k in range(len(u)) if u[k]], default=-1)
    if top < 0:
        return 'ALL'
    if top == 0:
        return []
    if top == 1:                       # u0 + u1*t
        return [-u[0]//u[1]] if u[0] % u[1] == 0 else []
    if top == 2:                       # u0 + u1*t + u2*t(t-1)/2  -> 2u0 + (2u1-u2)t + u2 t^2
        A, B, Cc = u[2], 2*u[1]-u[2], 2*u[0]
        D_ = B*B - 4*A*Cc
        if D_ < 0:
            return []
        import math
        r = math.isqrt(D_)
        if r*r != D_:
            return []
        out = []
        for s_ in (-B+r, -B-r):
            if s_ % (2*A) == 0:
                out.append(s_//(2*A))
        return sorted(set(out))
    return None

def exact_pins(CF, EX, D, log, w, v):
    """resolve every exact condition to a pin on t_w and/or t_v.
    returns (setw, setv, ok).  setw/setv are None when unconstrained."""
    setw = setv = None
    for a in EX:
        cf = CF[a]
        depw = any(cf[k][l] for k in range(1, D+1) for l in range(D+1))
        depv = any(cf[k][l] for k in range(D+1) for l in range(1, D+1))
        if not depw and not depv:
            if cf[0][0] != 0:
                return None, None, False          # constant nonzero: unfixable on this pair
            continue
        if depw and depv:
            log('        exact condition depends on BOTH wires -- not implemented, pair skipped')
            return None, None, False
        if depw:
            rs = newton_int_roots([cf[k][0] for k in range(D+1)])
            tgt = 'w'
        else:
            rs = newton_int_roots([cf[0][l] for l in range(D+1)])
            tgt = 'v'
        if rs is None:
            log('        exact condition of degree >2 -- not implemented, pair skipped')
            return None, None, False
        if rs == 'ALL':
            continue
        if not rs:
            return None, None, False              # no integer root: unsatisfiable on this pair
        if tgt == 'w':
            setw = set(rs) if setw is None else (setw & set(rs))
            if not setw:
                return None, None, False
        else:
            setv = set(rs) if setv is None else (setv & set(rs))
            if not setv:
                return None, None, False
    return setw, setv, True

def transpose_cf(cf):
    """swap the roles of t_w and t_v in a Newton coefficient table (cf[k][l], k over t_w)."""
    D = len(cf)-1
    return [[cf[j][i] for j in range(D+1)] for i in range(D+1)]

def mod_tv_sets(CF, MOD, tw, D):
    """with t_w pinned, the divisibility conditions are univariate in t_v: per prime power the
    intersection of the allowed residues.  Exhaustive."""
    PP = {}
    for a in MOD:
        for q, e in factor(abs(SL[a])//p).items():
            PP[q] = max(PP.get(q, 0), e)
    per = []
    for q in sorted(PP):
        e = PP[q]; m = q**e
        cand = None
        for a in MOD:
            fa = factor(abs(SL[a])//p)
            if q not in fa:
                continue
            ma = q**min(e, fa[q])
            rs = tv_roots(CF[a], tw, ma, q, min(e, fa[q]), D)
            if ma < m:
                rs = set(b for b in range(m) if b % ma in rs)
            cand = rs if cand is None else (cand & rs)
            if not cand:
                return None
        per.append((m, sorted(cand if cand is not None else range(m))))
    return per

def single_probe(vv, a, w, base, ZERO0, newatoms, log):
    """One-wire attempt on atom a, used to DISCOVER collateral when the atom has too few wires to
    form a pair (the |S|=32 handle-less atoms have exactly one).  Same two guards."""
    i = E.residx[a]
    ys = []
    for t in range(6):
        ow = vv[w]; vv[w] = ow + p*t; ys.append(E.run(vv)[i]); vv[w] = ow
    d = [ys[:]]
    for k in range(5):
        d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
    cf = [d[k][0] for k in range(6)]
    if is_exact(a):
        rs = newton_int_roots(cf)
        if rs is None or rs == 'ALL' or not rs:
            return False
        cands = rs[:4]
    else:
        if any(y % p for y in ys):
            return False
        cands = C.roots_c([x//p for x in cf], abs(SL[a])//p)[:40]
    for t in cands:
        snap = vv[:]
        vv[w] += p*t
        y = E.run(vv)[i]
        if (y == 0) if is_exact(a) else (y % abs(SL[a]) == 0):
            n = nzcount(vv)
            if n < base:
                log('     *** single wire x%d += p*%d  global %d -> %d' % (w, t, base, n))
                return True
            rn = E.run(vv)
            for i2 in ZERO0:
                if rn[i2]:
                    newatoms[E.res[i2]] += 1
        vv[:] = snap
    return False


FROZEN = set()

def forced_exact_pass(vv, hl, log):
    """T35 -- the last resort, and the only pass that is allowed to raise the global count.

    A handle-less atom has no cofactor to absorb anything, so its condition is `R == 0` over Z.
    On its admitted wires R(t) is LINEAR, so the root is UNIQUE: the shift is FORCED, not chosen,
    and refusing it because the global nonzero count does not immediately drop is refusing the only
    value the wire may take.  So: apply it, verify `R == 0` by direct recomputation, then FREEZE the
    wire (remove it from SHIFT, which every downstream pass filters on) so nothing can move it back,
    and let the ordinary guarded passes clear whatever it broke.  The global guard is untouched for
    every other kind of step; this one is justified by uniqueness, not by improvement."""
    for a in hl:
        i = E.residx[a]
        ws = sorted(set(x for x in vars_of(E.atoms[a]) if x in SHIFT) |
                    set(x for x in atomvalvars[a] if x in SHIFT))
        ws = [w for w in ws if influences(vv, a, w)]
        if len(ws) != 1:
            log('   forced-exact: %s has %d admitted wires -- not forced, skipped'
                % (a[:44], len(ws)))
            continue
        w = ws[0]
        ys = []
        for t in range(6):
            ow = vv[w]; vv[w] = ow + p*t; ys.append(E.run(vv)[i]); vv[w] = ow
        d = [ys[:]]
        for k in range(5):
            d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
        cf = [d[k][0] for k in range(6)]
        top = max([k for k in range(6) if cf[k]], default=0)
        rs = newton_int_roots(cf)
        if top != 1 or not isinstance(rs, list) or len(rs) != 1:
            log('   forced-exact: %s deg=%d roots=%s -- root not unique, skipped'
                % (a[:44], top, rs))
            continue
        t = rs[0]
        before = nzcount(vv)
        vv[w] += p*t
        if E.run(vv)[i] != 0:
            vv[w] -= p*t
            log('   forced-exact: direct recomputation REFUTED the fitted root -- skipped')
            continue
        after = nzcount(vv)
        SHIFT.discard(w); FROZEN.add(w)
        log('   *** FORCED exact step  x%d += p*%d  zeroes %s  (unique root, recomputed); '
            'global %d -> %d; wire FROZEN' % (w, t, a[:44], before, after))
        return True
    return False


def joint_pair(vv, a0, base, log):
    """clear atom a0 with a joint two-wire shift.  The group grows by folding in whatever the
    global guard says the candidate broke, and the wire set grows with it.  Mixed exact/mod."""
    def wires_of(a):
        ws = set(x for x in vars_of(E.atoms[a]) if x in SHIFT) | \
             set(x for x in atomvalvars[a] if x in SHIFT)
        return [x for x in sorted(ws) if influences(vv, a, x)]
    kind = 'EXACT (handle-less)' if is_exact(a0) else 'c=%d' % (abs(SL[a0])//p)
    GROUP = [a0]
    r0 = E.run(vv)
    ZERO0 = set(i for i, x in enumerate(r0) if x == 0)
    log('  atom %s  %s' % (a0[:56], kind))
    tried_pairs = set()
    for rd in range(ROUNDS):
        WSET = sorted(set(x for a in GROUP for x in wires_of(a)))
        pairs = [pr for pr in itertools.combinations(WSET, 2) if pr not in tried_pairs]
        log('     round %d: |group|=%d (%d exact) wires=%d new pairs=%d'
            % (rd, len(GROUP), sum(1 for a in GROUP if is_exact(a)), len(WSET), len(pairs)))
        newatoms = collections.Counter()
        for w in WSET[:8]:
            if single_probe(vv, a0, w, base, ZERO0, newatoms, log):
                return True
        for w, v in pairs[:PAIRCAP]:
            tried_pairs.add((w, v))
            CF = {}; okfit = True
            for a in GROUP:
                cf = C.fit2(vv, E.residx[a], w, v)
                if cf is None or not all(
                        C.probe2(vv, E.residx[a], w, tw, v, tv)//p == C.peval2_exact(cf, tw, tv)
                        for tw, tv in ((7, 11), (1234, 99), (55555, 4))):
                    okfit = False; break
                CF[a] = cf
            if not okfit:
                continue
            D = len(CF[a0])-1
            EX = [a for a in GROUP if is_exact(a)]
            MOD = [a for a in GROUP if not is_exact(a)]
            cands = []
            if EX:
                setw, setv, ok = exact_pins(CF, EX, D, log, w, v)
                if not ok:
                    continue
                if setw is not None and setv is not None:
                    cands = [(x, y) for x in sorted(setw)[:6] for y in sorted(setv)[:6]]
                elif setw is not None:
                    for x in sorted(setw)[:4]:
                        per = mod_tv_sets(CF, MOD, x, D) if MOD else []
                        if per is None:
                            continue
                        if not per:
                            cands.append((x, 0)); continue
                        mods = [m for m, _ in per]; sets = [rs for _, rs in per]
                        for _ in range(TRIES):
                            y = crt_list([(rs[rnd.randrange(len(rs))], m)
                                          for rs, m in zip(sets, mods)])
                            if y is not None:
                                cands.append((x, y))
                elif setv is not None:
                    # T35: the SYMMETRIC case.  This branch used to leave t_w at 0, so the
                    # divisibility conditions were never solved and the only candidate offered
                    # was the bare exact pin -- exactly the shift the global guard had already
                    # refused.  That is why |S|=32 stalled.  Transposing the coefficient table
                    # swaps the roles of the two wires and reuses the same exhaustive machinery.
                    CFT = {a: transpose_cf(CF[a]) for a in MOD}
                    for y in sorted(setv)[:4]:
                        per = mod_tv_sets(CFT, MOD, y, D) if MOD else []
                        if per is None:
                            continue
                        if not per:
                            cands.append((0, y)); continue
                        mods = [m for m, _ in per]; sets = [rs for _, rs in per]
                        for _ in range(TRIES):
                            x = crt_list([(rs[rnd.randrange(len(rs))], m)
                                          for rs, m in zip(sets, mods)])
                            if x is not None:
                                cands.append((x, y))
                else:
                    cands = []
            if not EX or not cands:
                if EX and not cands:
                    continue
                PP = {}
                for a in MOD:
                    for q, e in factor(abs(SL[a])//p).items():
                        PP[q] = max(PP.get(q, 0), e)
                per = []; dead = None
                for q in sorted(PP):
                    rs, exh, m = joint_rootsets(CF, MOD, q, PP[q])
                    per.append((m, rs))
                    if not rs:
                        dead = (m, exh); break
                if dead:
                    log('     (x%-6d,x%-6d) |group|=%d NO JOINT ROOT mod %d %s'
                        % (w, v, len(GROUP), dead[0], '(exhaustive)' if dead[1] else '(sampled)'))
                    continue
                mods = [m for m, _ in per]; sets = [rs for _, rs in per]
                seen = set()
                for _ in range(TRIES):
                    pick = [rs[rnd.randrange(len(rs))] for rs in sets]
                    tw = crt_list([(x, m) for (x, _), m in zip(pick, mods)])
                    tv = crt_list([(y, m) for (_, y), m in zip(pick, mods)])
                    if tw is None or tv is None or (tw, tv) in seen:
                        continue
                    seen.add((tw, tv)); cands.append((tw, tv))
            cleared = 0
            for tw, tv in cands[:TRIES]:
                if (tw, tv) == (0, 0):
                    continue
                good = True
                for a in GROUP:
                    y = C.probe2(vv, E.residx[a], w, tw, v, tv)
                    if (y != 0) if is_exact(a) else (y % abs(SL[a]) != 0):
                        good = False; break
                if not good:
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
            if cleared:
                log('     (x%-6d,x%-6d) |group|=%d %d candidate(s) verified, guard rejected all'
                    % (w, v, len(GROUP), cleared))
        add = [a for a, _ in newatoms.most_common() if a not in GROUP][:3]
        if not add:
            log('     no addable collateral -- giving up on %s' % a0[:44])
            break
        log('     folding in collateral: %s' % [a[:44] for a in add])
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

    def dump(why):
        """T38: dump after EVERY outer round, not only at the end.  close_T128.json was lost when
        that process was killed mid-round after 40 minutes -- a long run must always leave a
        scorable artefact, and a log is not an artefact."""
        json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
                  open(os.path.join(T, 'close_%s.json' % tag), 'w'))
        log('   [dumped close_%s.json after %s]' % (tag, why))

    for outer in range(outer_max):
        base = nzcount(vv); r = E.run(vv); gen += 1
        viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0]
        hl0 = [a for a in E.res if r[E.residx[a]] and a not in SL]
        log('outer %d: global nonzero %d, violated c-conditions %d, nonzero handle-less %d'
            % (outer, base, len(viol), len(hl0)))
        dump('outer %d (global nonzero %d)' % (outer, base))
        if not viol and not hl0:
            break
        if not viol:                       # only handle-less exact conditions left
            if handleless_pass(vv, base, log):
                continue
            if any(joint_pair(vv, a, base, log) for a in hl0):
                continue
            if forced_exact_pass(vv, hl0, log):
                continue
            log('   handle-less atoms remain and nothing moves them -> stop')
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
        hl = [a for a in E.res if r[E.residx[a]] and a not in SL]      # handle-less, exact-zero
        if not viol and not hl:
            log('   only the two TARGET congruences remain -- CLOSED'); break
        if not any(joint_pair(vv, a, base, log) for a in viol+hl):
            if hl and forced_exact_pass(vv, hl, log):
                gen += 1; continue
            log('   joint two-wire pass also stalled -> stop'); break
        gen += 1
    relift(vv); r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
              open(os.path.join(T, 'close_%s.json' % tag), 'w'))
    return nz

if __name__ == '__main__':
    tag = sys.argv[1]; spec = sys.argv[2]
    # 4th arg = ON-set SEED.  The historical convention is random.Random(7), whose samples at
    # n = 32 / 64 / 128 are NESTED PREFIXES of one chain -- ONE correlated sample, not three.
    # Any other seed draws an INDEPENDENT ON-set.
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    if ',' in spec:
        S = [int(x) for x in spec.split(',')]
    else:
        n = int(spec); r7 = random.Random(seed)
        S = [24601, 2081] if (n == 2 and seed == 7) else r7.sample(M['live'], n)
    omax = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    lf = open(os.path.join(T, 't_close2wj_%s.log' % tag), 'w')
    lf.write('S(seed=%d) = %s\n' % (seed, S)); print('S(seed=%d) = %s' % (seed, S), flush=True)
    t0 = time.time(); nz = close(S, tag, outer_max=omax, logf=lf); el = time.time()-t0
    msg = ('|S|=%-3d %-6s  NONZERO ATOMS = %d of %d   WALL = %.1f s  -> close_%s.json'
           % (len(S), tag, len(nz), len(E.res), el, tag))
    print(msg, flush=True); lf.write(msg+'\n')
    for a in nz[:10]:
        print('    ', a[:110], flush=True); lf.write('     '+a[:110]+'\n')
    lf.close()
