#!/usr/bin/env python3
"""agent V -- V8.  DOES THE HANDLE-LESS RESIDUE ACTUALLY BLOCK, OR IS IT A GUARD ARTIFACT?

v_hreach measured, at |S|=32, that each handle-less nonzero atom has an EXACT integer root on its
single influencing wire, and that taking the root breaks exactly ONE other atom -- and that atom is
a HANDLED one (c==1 in one case, c>1 in the other).  closeS4 refuses both because its guard is
"the global nonzero count must strictly decrease", and a 1-for-1 trade does not decrease it.

But a handled atom is not the same kind of object as a handle-less one:
  * c == 1  is absorbed by `relift` outright -- it is not a residue at all;
  * c >  1  is a divisibility condition, which is exactly what the existing solver discharges.
A handle-less atom has nothing to absorb into and MUST be exactly zero.

So the trade is not neutral, and the guard that scores it as neutral is measuring the wrong thing.
This file replaces it with a LEXICOGRAPHIC guard --

      (handle-less nonzero count, then global nonzero count)  must strictly decrease

-- and tests whether that closes |S|=32.  Every root is verified by direct recomputation, the
result is dumped, and it must be put in front of checker.py: nothing here is a claim until it is.

Usage: python3 v_hclose.py <tag> <n | comma,list>
"""
import os, sys, json, time, collections
from math import gcd
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, NV = B.E, B.SL, B.p, B.NV
atomh = B.GL['atomh']
HLS = set(a for a in E.res if len(atomh[a]) == 0)


def counts(vv):
    """(handle-less nonzero, global nonzero) AFTER relift -- relift absorbs what it can, so this
    measures what actually survives rather than what is transiently nonzero."""
    B.relift(vv)
    r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    return sum(1 for a in nz if a in HLS), len(nz), nz


def fit1(vv, a, w, D=4):
    i = E.residx[a]
    ys = [B.probe(vv, i, [w], [t]) for t in range(D+1)]
    d = [ys[:]]
    for k in range(D):
        d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
    cf = [d[k][0] for k in range(D+1)]
    return cf, max([k for k in range(D+1) if cf[k]], default=0)


def exact_roots(cf, top):
    if top == 0:
        return [0] if cf[0] == 0 else []
    if top == 1:
        return [-cf[0]//cf[1]] if cf[1] and cf[0] % cf[1] == 0 else []
    out = []
    for t in range(-8192, 8193):
        v, b = 0, 1
        for k, c in enumerate(cf):
            if k:
                b = b*(t-k+1)//k
            v += c*b
        if v == 0:
            out.append(t)
    return out


def handleless_pass(vv, verbose=True):
    """zero handle-less atoms by exact integer roots, under the LEXICOGRAPHIC guard."""
    prog = 0
    for _ in range(12):
        h0, g0, nz = counts(vv)
        hl = [a for a in nz if a in HLS]
        if not hl:
            break
        moved = False
        for a in hl:
            ws = sorted(w for w in B.wires_of(a) if B.influences(vv, a, w))
            for w in ws:
                cf, top = fit1(vv, a, w)
                for t in exact_roots(cf, top):
                    if B.probe(vv, E.residx[a], [w], [t]) != 0:
                        continue                      # GUARD 1: exact zero by recomputation
                    old = vv[w]
                    vv[w] = old + p*t
                    h1, g1, _ = counts(vv)
                    if (h1, g1) < (h0, g0):           # GUARD 2: lexicographic
                        if verbose:
                            print('   x%-6d += p*%-20d  (handle-less %d->%d, global %d->%d)  %s'
                                  % (w, t, h0, h1, g0, g1, a[:52]), flush=True)
                        prog += 1
                        moved = True
                        h0, g0 = h1, g1
                        break
                    vv[w] = old
                if moved:
                    break
            if moved:
                break
        if not moved:
            break
    return prog


def cond_pass(vv, verbose=True):
    """L's single-wire c>1 pass, with the same lexicographic guard so it cannot trade a divisibility
    condition for a handle-less nonzero."""
    gen = 0
    prog = 0
    for outer in range(10):
        h0, g0, _ = counts(vv)
        r = E.run(vv)
        gen += 1
        viol = B.violated(vv, r)
        if not viol:
            break
        wires = collections.defaultdict(list)
        for a in viol:
            for w in B.wires_of(a):
                wires[w].append(a)
        p0 = prog
        for w, ats in sorted(wires.items(), key=lambda kv: -len(kv[1])):
            Vt = [a for a in ats if B.influences(vv, a, w)]
            if not Vt:
                continue
            t = B.GL['solve_group3'](vv, Vt, w, gen, g0)
            if t:
                prog += 1
                h0, g0, _ = counts(vv)
                gen += 1
        if prog == p0:
            break
    return prog


if __name__ == '__main__':
    tag = sys.argv[1]
    arg = sys.argv[2]
    S = [int(x) for x in arg.split(',')] if ',' in arg else B.onset(int(arg))
    open(os.path.join(V, 'run_%s.pid' % tag), 'w').write('%d\n' % os.getpid())
    print('pid %d ; |S|=%d' % (os.getpid(), len(S)), flush=True)
    t0 = time.time()
    vv = B.greedy_init(S)
    h, g, _ = counts(vv)
    print('greedy fixpoint: handle-less nonzero %d, global nonzero %d' % (h, g), flush=True)
    for rd in range(8):
        a1 = cond_pass(vv)
        h, g, _ = counts(vv)
        print('round %d: c>1 pass accepted %d  -> handle-less %d, global %d' % (rd, a1, h, g),
              flush=True)
        a2 = handleless_pass(vv)
        h, g, _ = counts(vv)
        print('round %d: handle-less pass accepted %d  -> handle-less %d, global %d'
              % (rd, a2, h, g), flush=True)
        if a1 + a2 == 0:
            break
    h, g, nz = counts(vv)
    out = os.path.join(V, 'close_%s.json' % tag)
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]}, open(out, 'w'))
    print('\n|S|=%-4d %-8s NONZERO ATOMS = %d of 9032  (handle-less %d)   WALL = %.1f s'
          % (len(S), tag, g, h, time.time()-t0), flush=True)
    for a in nz:
        mark = 'HANDLE-LESS' if a in HLS else ('c=%d' % (abs(SL[a])//p) if SL.get(a) and SL[a] % p == 0 else '')
        print('     %-78s %s' % (a[:78], mark), flush=True)
    print('-> %s   ** run checker.py on it **' % out, flush=True)
