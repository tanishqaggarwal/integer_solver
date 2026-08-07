#!/usr/bin/env python3
"""agent V -- V7.  REACHABILITY WITH COLLATERAL PRICED.

v_hwho established that every handle-less atom seen nonzero has a wire on which R(t) is LINEAR with
an EXACT integer root, verified by recomputation.  That answers "is there a root".  It does not
answer "can the solver take it", which is a different question and the one T hit: an exact root
that breaks something else is refused by the global guard.

So price each root: apply it, recompute everything, and report

    * the global nonzero count before and after,
    * exactly which atoms the shift breaks and which it fixes,
    * whether the breakage is itself handle-less (fatal -- nothing absorbs it) or handled with
      c == 1 or c > 1 (absorbable / a divisibility condition the existing machinery can chase).

Then, for the atoms whose single root is refused, try the natural escalation: does moving a SECOND
wire jointly clear it?  Handle-less atoms need R == 0 exactly over Z, so this is an integer
system, not a congruence -- for a bilinear R the second wire's value is forced, and we test it
directly rather than searching.

Usage: python3 v_hreach.py <n>
"""
import os, sys, json, time, collections, itertools
from math import gcd
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, NV = B.E, B.SL, B.p, B.NV
atomh = B.GL['atomh']
HLS = set(a for a in E.res if len(atomh[a]) == 0)


def cls(a):
    """what kind of atom is this, for pricing breakage."""
    if a in HLS:
        return 'HANDLE-LESS (must be exactly 0)'
    s = SL.get(a)
    if not s:
        return 'handled, zero slope'
    c = abs(s)//p if s % p == 0 else None
    if c is None:
        return 'handled, slope not p-divisible'
    return 'handled c==1' if c == 1 else 'handled c>1 (c=%d)' % c


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


if __name__ == '__main__':
    n = int(sys.argv[1])
    tag = 'reach%d' % n
    open(os.path.join(V, 'run_%s.pid' % tag), 'w').write('%d\n' % os.getpid())
    S = B.onset(n)
    print('pid %d ; |S|=%d' % (os.getpid(), n), flush=True)
    t0 = time.time()
    vv = B.greedy_init(S)
    r0 = E.run(vv)
    base = [i for i, x in enumerate(r0) if x]
    hl = [E.res[i] for i in base if E.res[i] in HLS]
    print('greedy fixpoint: global nonzero %d, handle-less nonzero %d' % (len(base), len(hl)),
          flush=True)
    Z0 = set(i for i, x in enumerate(r0) if x == 0)

    out = {}
    for a in hl:
        print('\n=== %s ===' % a[:90], flush=True)
        ws = sorted(w for w in B.wires_of(a) if B.influences(vv, a, w))
        print('  influencing SHIFT wires: %d %s' % (len(ws), ['x%d' % w for w in ws]), flush=True)
        rec = {'wires': ws, 'roots': []}
        for w in ws:
            cf, top = fit1(vv, a, w)
            rs = exact_roots(cf, top)
            print('  x%-6d degree %d ; exact integer roots: %s' % (w, top, rs[:4]), flush=True)
            for t in rs:
                if B.probe(vv, E.residx[a], [w], [t]) != 0:
                    print('     t=%d REJECTED by recomputation (fit disagreed)' % t, flush=True)
                    continue
                old = vv[w]
                vv[w] = old + p*t
                r1 = E.run(vv)
                nz1 = [i for i, x in enumerate(r1) if x]
                broke = [E.res[i] for i in nz1 if i in Z0]
                fixed = [E.res[i] for i in base if not r1[i]]
                vv[w] = old
                bc = collections.Counter(cls(x) for x in broke)
                print('     t=%-16d global %d -> %d   fixed %d, broke %d  %s'
                      % (t, len(base), len(nz1), len(fixed), len(broke), dict(bc)), flush=True)
                for x in broke[:6]:
                    print('        breaks: %-64s [%s]' % (x[:64], cls(x)), flush=True)
                rec['roots'].append({'wire': w, 't': t, 'degree': top,
                                     'global_before': len(base), 'global_after': len(nz1),
                                     'broke': broke, 'fixed': fixed,
                                     'broke_classes': dict(bc),
                                     'accepted_by_global_guard': len(nz1) < len(base)})
        out[a] = rec

    json.dump({'S': S, 'global_nz': len(base), 'handleless_nz': hl, 'detail': out},
              open(os.path.join(V, 'v_hreach_%d.json' % n), 'w'), indent=1)
    print('\nwall %.1f s -> v_hreach_%d.json' % (time.time()-t0, n), flush=True)
