#!/usr/bin/env python3
"""perstep.py -- honest per-step / per-window breakdown, squaring-aware.

One comb "step" (group addition, windows 1..M-1) is exactly THREE field
multiplies plus glue:

    d  = x2 - x1            (linear word)
    e  = y2 - y1            (linear word)
    lam*d == e             GENERAL modmul          <- atom, general
    d != 0, d != p         two not_equal gadgets
    lam*lam == x3+x1+x2    SQUARING                 <- atom, squaring (cheaper)
    m  = x1 - x3           (linear word)
    lam*(x1-x3) == y3+y1   GENERAL modmul          <- atom, general

so per step = 2 general + 1 squaring, NOT 3 general.  A window additionally
carries two one-hot table look-ups (2^w entries each) and its one-hot selector
chain.  This isolates every piece at s=256 with the real p.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.normpath(os.path.join(HERE, '..', '..', 'squeeze'))
sys.path.insert(0, SQ)
import random                                                     # noqa: E402
from mmqb import MMQB, best_split, bin_split                      # noqa: E402
from mm import build_modmul                                       # noqa: E402
from ladder2 import lin_word, not_equal, _qword                   # noqa: E402

P = 2 ** 256 - 2 ** 32 - 977
S = 256


def measure(fn):
    Q = MMQB(chunk=16, mode='wallace')
    base = fn(Q)
    Q.finalize()
    st = Q.stats()
    return dict(vars=Q.n - base, clique=st['max_clique'],
                jbits=st['dynamic_range_bits'], couplers=st['couplers'])


def c_modmul(square):
    def f(Q):
        A = Q.mkword('A', S, lambda wv: 0)
        B = A if square else Q.mkword('B', S, lambda wv: 0)
        C = Q.mkword('C', S, lambda wv: 0)
        base = Q.n
        build_modmul(Q, P, A, B, C, mult='karatsuba', leaf=24, red='naf')
        return base
    return f


def c_linword(Q):
    x = Q.mkword('x', S, lambda wv: 0)
    y = Q.mkword('y', S, lambda wv: 0)
    base = Q.n
    lin_word(Q, P, 'd', [(1, x), (-1, y)], 0, lambda wv: 0, True)
    return base


def c_neq(Q):
    d = Q.mkword('d', S, lambda wv: 0)
    base = Q.n
    not_equal(Q, P, d, 0, 'ne0')
    not_equal(Q, P, d, P, 'nep')
    return base


def c_select(w):
    def f(Q):
        D = 1 << w
        rnd = random.Random(1)
        consts = [rnd.randrange(P) for _ in range(D)]
        sel = [Q.new(f"u{t}", 'input') for t in range(D)]
        prev = None
        for t, v in enumerate(sel[:-1]):
            pv = Q.new(f"pp{t}", 'onehot')
            lf = {pv: -1, v: 1}
            if prev is not None:
                lf[prev] = 1
            Q.add_square(lf, 0)
            prev = pv
        Q.add_square({prev: 1, sel[-1]: 1}, -1)
        base = Q.n
        W = Q.mkword('sx', S, lambda wv: 0)
        q, qlo = _qword(Q, P, 'sx', -((1 << S) - 1), max(consts), lambda wv: 0)
        tl = []
        for t in range(D):
            for sg, sh in best_split(consts[t]):
                tl.append(((sel[t],), sg, sh))
        for t, b in enumerate(W.bits):
            tl.append(((b,), -1, t))
        for t, b in enumerate(q.bits):
            for sg, sh in best_split(P):
                tl.append(((b,), -sg, sh + t))
        Q.assert_terms(tl, bin_split(-P * qlo), 'sel')
        return base
    return f


if __name__ == '__main__':
    gen = measure(c_modmul(False))
    sq = measure(c_modmul(True))
    lin = measure(c_linword)
    neq = measure(c_neq)
    print("PER-STEP COMPONENTS (s=256, real p, karatsuba(24)/naf/wallace)")
    print(f"{'component':>34} {'logical':>10} {'K':>3} {'|J|':>5} {'couplers':>12}")
    print("-" * 70)

    def show(name, r, mult=1):
        print(f"{name:>34} {r['vars']*mult:10,d} {r['clique']:3d} "
              f"2^{r['jbits']:<3d} {r['couplers']*mult:12,d}")
    show("GENERAL modmul  (lam*d)", gen)
    show("GENERAL modmul  (lam*(x1-x3))", gen)
    show("SQUARING        (lam*lam)", sq)
    show("linear word (d / e / m)  x3", lin, 3)
    show("not_equal x2 (d!=0, d!=p)", neq)
    step_vars = 2 * gen['vars'] + sq['vars'] + 3 * lin['vars'] + neq['vars']
    step_coup = 2 * gen['couplers'] + sq['couplers'] + 3 * lin['couplers'] + neq['couplers']
    print("-" * 70)
    print(f"{'STEP (2 general + 1 squaring + glue)':>34} {step_vars:10,d} {'5':>3} "
          f"{'2^9':>5} {step_coup:12,d}")
    print(f"\n  multiply mix per step: 2 GENERAL ({gen['vars']:,d} each) + "
          f"1 SQUARING ({sq['vars']:,d})")
    print(f"  if all three were counted as general: "
          f"{3*gen['vars']:,d} vars in multiplies")
    print(f"  squaring-aware actual multiplies:     "
          f"{2*gen['vars']+sq['vars']:,d} vars  "
          f"(saves {gen['vars']-sq['vars']:,d} = "
          f"{100*(gen['vars']-sq['vars'])/(3*gen['vars']):.1f}% of multiply cost)")
    print("\nPER-WINDOW = step + 2 table look-ups + selector chain (grows with 2^w):")
    print(f"{'w':>4} {'2 selects+chain':>18} {'window total':>14}")
    for w in (4, 6, 8, 10, 12):
        selr = measure(c_select(w))
        win = step_vars + 2 * selr['vars']
        print(f"{w:4d} {2*selr['vars']:18,d} {win:14,d}")
