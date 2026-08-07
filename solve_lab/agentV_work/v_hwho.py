#!/usr/bin/env python3
"""agent V -- V5.  WHICH handle-less atoms can go nonzero, and are they reachable?

Only 33 of the 5,351 handle-less atoms have ANY SHIFT wire in scope; the other 5,318 cannot be
moved by any shift the solver has.  So the interesting question splits cleanly:

  * an atom that is never nonzero is not a gap however unreachable it is;
  * an atom that goes nonzero AND has no shift wire is a genuine dead end for THIS knob set;
  * an atom that goes nonzero and has a wire is a solver-coverage question -- does R(t) have an
    exact integer root on that wire?  ("exact" = R(t) == 0 over Z, not R(t) == 0 mod anything,
    because a handle-less atom has nothing to absorb a remainder into.)

Usage: python3 v_hwho.py
"""
import os, sys, json, collections, time, random
from math import gcd
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, SHIFT, NV = B.E, B.SL, B.p, B.SHIFT, B.NV
g = B.GL
atomh = g['atomh']; M = g['M']; vars_of = g['vars_of']; fa = g['fa']; handle = g['handle']

HL = [a for a in E.res if len(atomh[a]) == 0]
HLS = set(HL)


def exact_root_on(vv, a, w, D=4):
    """is there an integer t with R_a(t) EXACTLY 0 when x_w += p*t?  Fit, solve, then VERIFY by
    direct recomputation.  Returns t or None."""
    i = E.residx[a]
    ys = [B.probe(vv, i, [w], [t]) for t in range(D+1)]
    d = [ys[:]]
    for k in range(D):
        d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
    cf = [d[k][0] for k in range(D+1)]
    top = max([k for k in range(D+1) if cf[k]], default=0)
    cands = []
    if top == 0:
        return (0, top) if cf[0] == 0 else None
    if top == 1:
        # R(t) = cf0 + cf1*t  -> exact integer root iff cf1 | cf0
        if cf[1] and cf[0] % cf[1] == 0:
            cands.append(-cf[0]//cf[1])
    else:
        # degree >= 2: any integer root divides the constant term (rational root theorem on the
        # monomial form); enumerate small |t| and the divisors of R(0)
        for t in range(-4096, 4097):
            if _neval(cf, t) == 0:
                cands.append(t)
    for t in cands:
        if B.probe(vv, i, [w], [t]) == 0:          # GUARD: direct recomputation, exact zero
            return (t, top)
    return None


def _neval(cf, t):
    tot = 0; b = 1
    for k, c in enumerate(cf):
        if k:
            b = b*(t-k+1)//k
        tot += c*b
    return tot


if __name__ == '__main__':
    t0 = time.time()
    withshift = {a: sorted(B.wires_of(a)) for a in HL if B.wires_of(a)}
    print('handle-less atoms with a SHIFT wire in scope: %d of %d' % (len(withshift), len(HL)),
          flush=True)
    for a, ws in sorted(withshift.items(), key=lambda kv: -len(kv[1])):
        print('   %-64s wires %s' % (a[:64], ['x%d' % w for w in ws]), flush=True)

    print('\n--- which handle-less atoms actually go nonzero, over L\'s Random(7) ON-sets ---',
          flush=True)
    hits = collections.defaultdict(list)
    rows = []
    for n in [1, 2, 4, 8, 16, 17, 24, 32, 48, 64]:
        S = B.onset(n)
        vv = B.greedy_init(S)
        r = E.run(vv)
        nz = [E.res[i] for i, x in enumerate(r) if x]
        hl = [a for a in nz if a in HLS]
        rows.append((n, len(nz), len(hl)))
        for a in hl:
            hits[a].append(n)
        print('   |S|=%-4d global nonzero %-4d handle-less nonzero %-3d %s'
              % (n, len(nz), len(hl), [x[:44] for x in hl]), flush=True)

    print('\n--- reachability of the ones that DO go nonzero ---', flush=True)
    verdicts = {}
    for a, sizes in sorted(hits.items()):
        ws = sorted(B.wires_of(a))
        print('\n  %s' % a[:100], flush=True)
        print('     nonzero at |S| = %s ; SHIFT wires in scope: %d %s'
              % (sizes, len(ws), ['x%d' % w for w in ws][:8]), flush=True)
        if not ws:
            print('     >>> NO SHIFT WIRE AT ALL -- unreachable for this knob set', flush=True)
            verdicts[a] = {'wires': 0, 'reachable': False, 'why': 'no shift wire in scope'}
            continue
        S = B.onset(sizes[-1])
        vv = B.greedy_init(S)
        found = []
        for w in ws:
            if not B.influences(vv, a, w):
                continue
            res = exact_root_on(vv, a, w)
            if res:
                found.append((w, res[0], res[1]))
                print('     x%-6d : EXACT integer root t=%d (degree %d), verified by '
                      'recomputation' % (w, res[0], res[1]), flush=True)
            else:
                print('     x%-6d : no exact integer root' % w, flush=True)
        verdicts[a] = {'wires': len(ws), 'reachable': bool(found), 'roots': found}

    json.dump({'with_shift_wire': {a: ws for a, ws in withshift.items()},
               'rows': rows,
               'hits': {a: v for a, v in hits.items()},
               'verdicts': verdicts},
              open(os.path.join(V, 'v_hwho.json'), 'w'), indent=1)
    print('\nwall %.1f s -> v_hwho.json' % (time.time()-t0), flush=True)
