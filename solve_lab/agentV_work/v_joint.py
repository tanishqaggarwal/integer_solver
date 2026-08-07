#!/usr/bin/env python3
"""agent V -- the GENERAL k-wire joint solve.

L's `closeS4` moves ONE wire at a time and accepts a shift only if the global nonzero-atom count
strictly decreases.  A residue that needs several wires moved TOGETHER is indistinguishable, under
that search, from a residue with no solution at all.  This module removes that restriction for
arbitrary k.

Method, per component (a set of coupled conditions and the wire set they live on):
  1. fit each atom's residual/p as an EXACT k-variate Newton polynomial (integer, no division);
     VALIDATE the fit against direct recomputation at random points before using it;
  2. for each prime power q^e dividing any modulus in the component: enumerate the shift on k-1
     wires and ROOT-FIND the last from the resulting univariate polynomial -- never a k-fold loop;
  3. CRT the per-prime-power solution sets together;
  4. accept only after (a) DIRECT RECOMPUTATION of every atom in the group and (b) the global
     nonzero-atom count strictly decreasing.

Cost is q^(e(k-1)) per prime power, so the affordable k is set by the LARGEST prime power in the
component's moduli, not by k alone.  v_cost.py prices that; this file measures it.

Usage: python3 v_joint.py <tag> <n | comma,list> [kmax] [budget]
"""
import os, sys, json, time, itertools, collections, random
from math import gcd
import v_base as B
import v_poly as P

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, SHIFT, NV = B.E, B.SL, B.p, B.SHIFT, B.NV

DEG = 4                  # fit degree per variable; L's bound is 3, we fit one higher and CHECK
PPCAP = 4_000_000        # per-prime-power enumeration budget for k=2 (q^e above this is skipped)
rnd = random.Random(20260807)

STATS = collections.Counter()


# ------------------------------------------------------------------ probing / fitting
def probe_k(vv, i, wires, ts):
    old = [vv[w] for w in wires]
    for w, t in zip(wires, ts):
        vv[w] += p*t
    y = E.run(vv)[i]
    for w, o in zip(wires, old):
        vv[w] = o
    return y


def fit_k(vv, a, wires, D=DEG):
    """exact k-variate Newton fit of residual(a)/p.  None if the grid is not p-divisible."""
    i = E.residx[a]
    bad = [False]

    def f(idx):
        y = probe_k(vv, i, wires, idx)
        if y % p:
            bad[0] = True
            return 0
        return y//p
    cf = P.nd_fit(f, len(wires), D)
    if bad[0]:
        return None
    return cf


def validate_fit(vv, a, wires, cf, npts=4):
    """the guard the lab insists on: never trust an expansion without direct recomputation."""
    i = E.residx[a]
    for _ in range(npts):
        ts = tuple(rnd.randrange(0, 10**5) for _ in wires)
        y = probe_k(vv, i, wires, ts)
        if y % p or y//p != P.nd_eval(cf, ts):
            return False
    return True


# ------------------------------------------------------------------ joint root sets mod q^e
def _reduce_last(cf, pre, D, m):
    """substitute the first k-1 shifts, leaving Newton coefficients in the last variable mod m."""
    Bp = []
    for t in pre:
        b = [1]
        for i in range(1, D+1):
            b.append(b[-1]*(t-i+1)//i)
        Bp.append([x % m for x in b])
    u = [0]*(D+1)
    for idx, c in cf.items():
        if not c:
            continue
        v = c % m
        for j in range(len(pre)):
            v = v*Bp[j][idx[j]] % m
            if not v:
                break
        if v:
            u[idx[-1]] = (u[idx[-1]] + v) % m
    return u


def joint_roots_pp(CF, MODS, wires, q, e, cap):
    """EXHAUSTIVE over the first k-1 wires mod q^e; root-find the last.  [] if none, None if the
    enumeration exceeds `cap` (reported, never silently skipped)."""
    k = len(wires)
    m = q**e
    n_outer = m**(k-1)
    if n_outer > cap:
        return None
    need = [a for a in CF if MODS[a] % q == 0]
    if not need:
        return 'ANY'
    ma = {a: q**min(e, B.factor(MODS[a])[q]) for a in need}
    D = max(max(idx) for idx in CF[need[0]])
    out = []
    for pre in itertools.product(range(m), repeat=k-1):
        cand = None
        for a in need:
            mm = ma[a]
            u = _reduce_last(CF[a], pre, D, mm)
            rs = P.roots_pp(u, q, B.factor(MODS[a])[q] if mm < m else e)
            rs = set(rs)
            if mm < m:
                rs = set(b for b in range(m) if b % mm in rs)
            cand = rs if cand is None else (cand & rs)
            if not cand:
                break
        if cand:
            out.extend(tuple(pre)+(b,) for b in sorted(cand))
    return out


# ------------------------------------------------------------------ the component solve
def solve_component(vv, target, wires, base, budget=400, cap=PPCAP, verbose=True):
    """try to clear `target` by moving `wires` TOGETHER, preserving every other c>1 condition any
    of those wires influences.  Returns the accepted shift tuple, or None."""
    k = len(wires)
    ws = list(wires)
    r = E.run(vv)
    coll = [a for a in B.CGT2 if a != target and (B.wires_of(a) & set(ws))
            and r[E.residx[a]] % abs(SL[a]) == 0
            and any(B.influences(vv, a, w) for w in ws)]
    GROUP = [target] + coll
    MODS = {a: abs(SL[a])//p for a in GROUP}
    if verbose:
        print('    group: target c=%d + %d collateral c>1 condition(s); wires %s'
              % (MODS[target], len(coll), ['x%d' % w for w in ws]), flush=True)
    CF = {}
    for a in GROUP:
        cf = fit_k(vv, a, ws, DEG)
        if cf is None:
            if verbose:
                print('    fit failed (not p-divisible on the grid) for %s' % a[:50], flush=True)
            return None
        if not validate_fit(vv, a, ws, cf):
            if verbose:
                print('    FIT MISMATCH vs recomputation for %s -- refusing' % a[:50], flush=True)
            return None
        degs = P.nd_degrees(cf)
        STATS['maxdeg_%d' % max(degs)] += 1
        if max(degs) > 3:
            print('    !! DEGREE %d observed (> L\'s bound of 3) on %s' % (max(degs), a[:50]), flush=True)
        CF[a] = cf

    PP = {}
    for a in GROUP:
        for q, e in B.factor(MODS[a]).items():
            PP[q] = max(PP.get(q, 0), e)
    per = []
    for q in sorted(PP):
        e = PP[q]
        t0 = time.time()
        rs = joint_roots_pp(CF, MODS, ws, q, e, cap)
        if rs is None:
            if verbose:
                print('    mod %-10d : SKIPPED, enumeration %d^%d exceeds cap %d'
                      % (q**e, q**e, k-1, cap), flush=True)
            return 'CAPPED'
        if rs == 'ANY':
            continue
        if verbose:
            print('    mod %-10d : %d joint %d-tuple(s)  [%.1fs, exhaustive over %d wire(s)]'
                  % (q**e, len(rs), k, time.time()-t0, k-1), flush=True)
        if not rs:
            if verbose:
                print('    >>> NO JOINT ROOT mod %d -- the group is simultaneously unsatisfiable '
                      'on these %d wires, exhaustively.' % (q**e, k), flush=True)
            return None
        per.append((q**e, rs))
    if not per:
        return None

    mods = [m for m, _ in per]
    sets = [rs for _, rs in per]
    total = 1
    for rs in sets:
        total *= len(rs)
    if verbose:
        print('    joint solution count mod lcm = %d ; sampling up to %d' % (total, budget), flush=True)
    seen = set()
    tried = cleared = 0
    order = list(itertools.product(*sets)) if total <= budget else None
    while tried < budget:
        if order is not None:
            if tried >= len(order):
                break                                 # exhaustive: the whole solution set is done
            pick = order[tried]
        else:
            pick = tuple(rs[rnd.randrange(len(rs))] for rs in sets)
        tried += 1
        ts = []
        okc = True
        for j in range(k):
            t = B.crt_list([(pick[s][j], mods[s]) for s in range(len(mods))])
            if t is None:
                okc = False
                break
            ts.append(t)
        if not okc or tuple(ts) in seen or all(t == 0 for t in ts):
            continue
        seen.add(tuple(ts))
        if not all(probe_k(vv, E.residx[a], ws, ts) % abs(SL[a]) == 0 for a in GROUP):
            continue                                  # GUARD 1: direct recomputation
        cleared += 1
        old = [vv[w] for w in ws]
        for w, t in zip(ws, ts):
            vv[w] += p*t
        n = B.nzcount(vv)
        if n < base:                                  # GUARD 2: global nonzero count
            if verbose:
                print('    *** ACCEPTED %s -> global %d -> %d  (%d sampled, %d cleared)'
                      % (', '.join('x%d += p*%d' % (w, t) for w, t in zip(ws, ts)),
                         base, n, tried, cleared), flush=True)
            return ts
        for w, o in zip(ws, old):
            vv[w] = o
    if verbose:
        print('    %d candidate(s) sampled, %d cleared the group by recomputation, '
              'none passed the global guard' % (tried, cleared), flush=True)
    return None


# ------------------------------------------------------------------ driver
def close(S, tag, kmax=3, budget=400, cap=PPCAP, outer_max=14):
    t0 = time.time()
    vv = B.greedy_init(S)
    print('greedy init %.1f s -> global nonzero %d' % (time.time()-t0, B.nzcount(vv)), flush=True)
    gen = 0
    for outer in range(outer_max):
        base = B.nzcount(vv)
        r = E.run(vv)
        gen += 1
        viol = B.violated(vv, r)
        print('outer %d: global nonzero %d, violated c-conditions %d' % (outer, base, len(viol)),
              flush=True)
        if not viol:
            break
        # ---- pass 1: single wire (cheap; L's solver, unchanged)
        wires = collections.defaultdict(list)
        for a in viol:
            for w in B.wires_of(a):
                wires[w].append(a)
        prog = 0
        for w, ats in sorted(wires.items(), key=lambda kv: -len(kv[1])):
            Vt = [a for a in ats if B.influences(vv, a, w)]
            if not Vt:
                continue
            t = B.GL['solve_group3'](vv, Vt, w, gen, base)
            if t:
                prog += 1
                base = B.nzcount(vv)
                gen += 1
        print('   single-wire pass: %d accepted' % prog, flush=True)
        if prog:
            continue
        # ---- pass 2: k wires jointly, k = 2 .. kmax
        r = E.run(vv)
        viol = [a for a in B.violated(vv, r) if not any(t in a for t in B.TGTW)]
        if not viol:
            print('   only target congruences remain -> STOP', flush=True)
            break
        for a in viol:
            ws = sorted(w for w in B.wires_of(a) if B.influences(vv, a, w))
            c = abs(SL[a])//p
            print('  condition c=%d = %s  on %d wire(s)'
                  % (c, '*'.join('%d^%d' % qe for qe in sorted(B.factor(c).items())), len(ws)),
                  flush=True)
            done = False
            for k in range(2, min(kmax, len(ws))+1):
                for comb in itertools.combinations(ws, k):
                    res = solve_component(vv, a, comb, base, budget, cap)
                    if res == 'CAPPED':
                        continue
                    if res:
                        prog += 1
                        base = B.nzcount(vv)
                        done = True
                        break
                if done:
                    break
            if done:
                break
        if prog == 0:
            print('   joint pass to k=%d also stalled -> STOP' % kmax, flush=True)
            break
    B.relift(vv)
    r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    out = os.path.join(V, 'close_%s.json' % tag)
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]}, open(out, 'w'))
    el = time.time()-t0
    print('\n|S|=%-4d %-8s NONZERO ATOMS = %d of 9032   WALL = %.1f s  -> %s'
          % (len(S), tag, len(nz), el, out), flush=True)
    for a in nz:
        print('     %s' % a[:110], flush=True)
    print('degree histogram over fitted atoms: %s'
          % {k_: v for k_, v in sorted(STATS.items())}, flush=True)
    return nz, vv


if __name__ == '__main__':
    tag = sys.argv[1]
    arg = sys.argv[2]
    S = [int(x) for x in arg.split(',')] if ',' in arg else B.onset(int(arg))
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    budget = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    # The job writes its OWN pid.  $! is unreliable through setsid/nohup wrappers (it names the
    # wrapper, which then exits), and matching a command line has failed four times in this lab.
    # This is the only identification that cannot name the wrong process.
    open(os.path.join(V, 'run_%s.pid' % tag), 'w').write('%d\n' % os.getpid())
    print('pid %d -> run_%s.pid' % (os.getpid(), tag), flush=True)
    print('S (|S|=%d) = %s' % (len(S), S), flush=True)
    close(S, tag, kmax, budget)
