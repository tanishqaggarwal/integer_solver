#!/usr/bin/env python3
"""window.py -- marginal cost of ONE comb window at s = 256, old vs new.

A window is: two one-hot table look-ups, four linear words, three modular
multiplications (one of them a squaring) and two "d != 0" gadgets.  x1, y1 are
carried in from the previous window and are subtracted out, exactly as
../resources.py:marginal_window does, so the numbers compose:

    ladder = ceil(256/w) * window(w)
"""
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mmqb import MMQB, best_split, bin_split                     # noqa: E402
from mm import build_modmul                                      # noqa: E402
from ladder2 import lin_word, not_equal, _qword                  # noqa: E402

P = 2 ** 256 - 2 ** 32 - 977
HW = 4400
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'windows.json')


def window(w, mode='wallace', mult='karatsuba', leaf=24, red='naf',
           naf_consts=True, chunk=16, p=P, seed=1, naf_merge=True):
    s = p.bit_length()
    D = 1 << w
    rnd = random.Random(seed)
    consts = [rnd.randrange(p) for _ in range(2 * D)]
    Q = MMQB(chunk=chunk, mode=mode, naf_merge=naf_merge)
    split = best_split if naf_consts else bin_split

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

    def table(name, off):
        W = Q.mkword(name, s, lambda wv: 0)
        q, qlo = _qword(Q, p, name, -((1 << s) - 1), max(consts), lambda wv: 0)
        tl = []
        for t in range(D):
            for sg, sh in split(consts[off + t]):
                tl.append(((sel[t],), sg, sh))
        for t, b in enumerate(W.bits):
            tl.append(((b,), -1, t))
        for t, b in enumerate(q.bits):
            for sg, sh in split(p):
                tl.append(((b,), -sg, sh + t))
        Q.assert_terms(tl, bin_split(-p * qlo), f"sel:{name}")
        return W

    t0 = time.time()
    x2 = table("x2", 0)
    y2 = table("y2", D)
    base_words = Q.n
    x1 = Q.mkword("x1", s, lambda wv: 0)
    y1 = Q.mkword("y1", s, lambda wv: 0)
    d = lin_word(Q, p, "d", [(1, x2), (-1, x1)], 0, lambda wv: 0, naf_consts)
    e = lin_word(Q, p, "e", [(1, y2), (-1, y1)], 0, lambda wv: 0, naf_consts)
    lam = Q.mkword("lam", s, lambda wv: 0)
    build_modmul(Q, p, lam, d, [(e, 1)], mult=mult, leaf=leaf, red=red, tag="lam")
    for c in (0, p):
        not_equal(Q, p, d, c, f"ne{c}")
    x3 = Q.mkword("x3", s, lambda wv: 0)
    build_modmul(Q, p, lam, lam, [(x3, 1), (x1, 1), (x2, 1)], mult=mult, leaf=leaf,
                 red=red, tag="x3e")
    mm = lin_word(Q, p, "m", [(1, x1), (-1, x3)], 0, lambda wv: 0, naf_consts)
    y3 = Q.mkword("y3", s, lambda wv: 0)
    build_modmul(Q, p, lam, mm, [(y3, 1), (y1, 1)], mult=mult, leaf=leaf, red=red,
                 tag="y3e")
    Q.finalize()
    st = Q.stats()
    return dict(vars=Q.n - 2 * s, clique=st['max_clique'],
                jbits=st['dynamic_range_bits'], couplers=st['couplers'],
                ands=st['and_vars'], carries=st['carry_bits'],
                secs=round(time.time() - t0, 1))


def phys(v, c):
    return v * max(1.0, c / 6.0)


if __name__ == '__main__':
    import math
    res = json.load(open(RES)) if os.path.exists(RES) else {}
    cfgs = []
    for w in (6, 7, 8, 9, 10):
        cfgs.append((f"OLD-STYLE school/quotient/wallace w={w}",
                     dict(w=w, mult='schoolbook', red='quotient', mode='wallace',
                          naf_consts=False, naf_merge=False)))
        cfgs.append((f"SQUEEZED karatsuba(24)/naf/wallace w={w}",
                     dict(w=w, mult='karatsuba', leaf=24, red='naf',
                          mode='wallace', naf_consts=True)))
    print(f"{'configuration':>46} {'window':>11} {'K':>4} {'|J|':>6} "
          f"{'windows':>8} {'LADDER':>14} {'physical':>15} {'vs 4400':>9}")
    print("-" * 122)
    for key, kw in cfgs:
        if key not in res:
            res[key] = window(**kw)
            json.dump(res, open(RES, 'w'), indent=1, sort_keys=True)
        r = res[key]
        M = math.ceil(256 / kw['w'])
        tot = M * r['vars']
        ph = phys(tot, r['clique'])
        print(f"{key:>46} {r['vars']:11,d} {r['clique']:4d} 2^{r['jbits']:<4d} "
              f"{M:8d} {tot:14,d} {ph:15,.0f} {ph/HW:8,.0f}x")
