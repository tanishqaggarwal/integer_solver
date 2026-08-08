#!/usr/bin/env python3
"""ATTACK 1c: isolated ground-state enumeration of the two ladder-only gadgets
the task flags: the one-hot sequential counter and the not_equal (d != c) gadget.
Also probes the {0,p} completeness of not_equal used for the d!=0 loophole.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.join(HERE, '..', '..', 'squeeze')
sys.path.insert(0, os.path.abspath(SQ))
from mmqb import MMQB
import ladder2
import verify


def test_onehot(D, mode='wallace'):
    """build ONLY the sequential-counter one-hot on D selectors, enumerate all E=0."""
    Q = MMQB(chunk=16, mode=mode)
    u = [Q.new(f"u_{t}", 'input') for t in range(D)]
    for t, v in enumerate(u):
        Q.trace.append(('word', f"u_{t}", [v], (lambda wv, t=t: 1 if wv['_u'] == t else 0)))
    prev = None
    for t, v in enumerate(u[:-1]):
        pv = Q.new(f"p_{t}", 'onehot')
        Q.trace.append(('word', f"p_{t}", [pv], (lambda wv, t=t: 1 if wv['_u'] <= t else 0)))
        lin = {pv: -1, v: 1}
        if prev is not None: lin[prev] = 1
        Q.add_square(lin, 0)
        prev = pv
    Q.add_square({prev: 1, u[-1]: 1}, -1)
    Q.finalize()
    zeros = verify.zero_states(Q)
    onehots = set()
    bad = []
    for x in zeros:
        assert Q.energy(x) == 0
        sel = tuple(x[u[t]] for t in range(D))
        onehots.add(sel)
        if sum(sel) != 1:
            bad.append(sel)
    truth = set(tuple(1 if i == t else 0 for i in range(D)) for t in range(D))
    ok = onehots == truth
    print(f"  one-hot D={D} {mode:7s}: |E=0 selector patterns|={len(onehots)} truth={len(truth)} "
          f"{'OK' if ok and not bad else '*** BAD: '+str(sorted(onehots-truth))+' ***'}")
    return ok and not bad


def test_notequal(s, c, mode='wallace'):
    """build ONLY not_equal(W, c) on an s-bit free word W; E=0 states must be exactly {W != c}."""
    p = (1 << s) - 1  # dummy modulus, unused by not_equal itself
    Q = MMQB(chunk=16, mode=mode)
    W = Q.mkword('W', s, lambda wv: wv['_w'])
    ladder2.not_equal(Q, p, W, c, 'ne')
    Q.finalize()
    zeros = verify.zero_states(Q)
    wvals = set()
    for x in zeros:
        assert Q.energy(x) == 0
        wvals.add(sum(x[W.bits[t]] << t for t in range(s)))
    truth = set(w for w in range(1 << s) if w != c)
    ok = wvals == truth
    print(f"  not_equal s={s} c={c} {mode:7s}: |E=0 W|={len(wvals)} truth={len(truth)} "
          f"c in E=0? {c in wvals}  {'OK' if ok else '*** BAD ***'}")
    if not ok:
        print("     spurious(E=0 but ==c):", sorted(wvals & {c}), " missing:", sorted(truth - wvals)[:5])
    return ok


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    allok = True
    print("=== one-hot sequential counter ===")
    for D in (2, 3, 4, 5):
        for mode in ('wallace', 'binary'):
            allok &= test_onehot(D, mode)
    print("\n=== not_equal gadget (probe c=0 and c=p analogues, all c) ===")
    for s in (3, 4, 5):
        for c in range(1 << s):
            allok &= test_notequal(s, c, 'wallace')
    print(f"\nALL GADGETS FAITHFUL: {allok}")
