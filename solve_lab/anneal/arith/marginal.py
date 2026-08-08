#!/usr/bin/env python3
"""marginal.py -- exact marginal cost of ONE comb window / ONE S_3 step at any s.

Same method and the same accounting convention as the existing resources.py:
build a single window against fresh x1,y1 accumulator registers and subtract the
2*s bits of those registers, which belong to the previous window.  Random p and
random table constants: the cost of this encoding depends only on the SIZES,
which is why resources.py does the same.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ladder import not_equal              # noqa: E402
from enc import Ladder2                   # noqa: E402


def _prime(s, rnd):
    return (1 << (s - 1)) + 2 * rnd.randrange(1 << (s - 3)) + 1


def window(s, w, mode='binary', mux=True, kdepth=0, kmin=8, signed=False,
           onehot='square', chunk=16, seed=1, toom=0):
    """one comb window: table look-ups + 3 modular multiplications + linear words."""
    rnd = random.Random(seed)
    p = _prime(s, rnd)
    D = (1 << (w - 1)) if signed else (1 << w)
    L = Ladder2(p, chunk=chunk, mode=mode, kdepth=kdepth, kmin=kmin, toom=toom)
    Q = L.qb
    zero = lambda wv: 0

    def word(name):
        return Q.word(name, s, zero)

    dw = (w - 1) if signed else w
    if onehot == 'tree':
        sel = L.onehot_tree("u", dw, "_u")
    else:
        sel = L.onehot_square("u", dw, "_u")
    cx = [rnd.randrange(p) for _ in range(D)]
    cy = [rnd.randrange(p) for _ in range(D)]
    names = [f"u{t}" for t in range(D)]
    x2 = L.sel_word("x2", sel, names, cx, zero, mux)
    ya = L.sel_word("y2a", sel, names, cy, zero, mux)
    if signed:
        sg = Q.new("sg", 'input')
        z = [Q.AND(sg, v) for v in ya]
        y2 = Q.word("y2", s, zero)
        poly = {}
        for t, v in enumerate(ya):
            poly[(v,)] = poly.get((v,), 0) + (1 << t)
        for t, v in enumerate(z):
            poly[(v,)] = poly.get((v,), 0) - 2 * (1 << t)
        for t, v in enumerate(y2):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        L.congruent(poly, 0, "lin:y2", zero)
    else:
        y2 = ya

    x1, y1 = word("x1"), word("y1")
    d = L.lin_word("d", [(1, x2, "x2"), (-1, x1, "x1")], 0, zero)
    e = L.lin_word("e", [(1, y2, "y2"), (-1, y1, "y1")], 0, zero)
    lam = word("lam")
    L.mul_eq("lam", lam, d, "lam", "d", [(1, e, "e")], 0)
    for c in (0, p):
        not_equal(L, Q, d, "d", c, f"ne{c}")
    x3 = word("x3")
    L.mul_eq("x3", lam, lam, "lam", "lam",
             [(1, x3, "x3"), (1, x1, "x1"), (1, x2, "x2")], 0)
    mm = L.lin_word("m", [(1, x1, "x1"), (-1, x3, "x3")], 0, zero)
    y3 = word("y3")
    L.mul_eq("y3", lam, mm, "lam", "m", [(1, y3, "y3"), (1, y1, "y1")], 0)
    Q.finalize()
    st = Q.stats()
    st['vars'] -= 2 * s          # x1, y1 belong to the previous window
    return st


def semaev_step(s, w, mode='binary', mux=True, kdepth=0, kmin=8,
                onehot='square', chunk=16, seed=1):
    """one S_3 chain step: one x-only table look-up + 6 modular multiplications."""
    rnd = random.Random(seed)
    p = _prime(s, rnd)
    Bc = rnd.randrange(p)
    D = 1 << (w - 1)
    L = Ladder2(p, chunk=chunk, mode=mode, kdepth=kdepth, kmin=kmin)
    Q = L.qb
    zero = lambda wv: 0

    def word(name):
        return Q.word(name, s, zero)

    if onehot == 'tree':
        sel = L.onehot_tree("u", w - 1, "_u")
    else:
        sel = L.onehot_square("u", w - 1, "_u")
    cx = [rnd.randrange(p) for _ in range(D)]
    xj = L.sel_word("xt", sel, [f"u{t}" for t in range(D)], cx, zero, mux)
    z = word("z")               # incoming accumulator x
    zn = word("zn")             # outgoing accumulator x
    v = L.lin_word("v", [(1, z, "z"), (1, xj, "xt")], 0, zero)
    dd = L.lin_word("dd", [(1, z, "z"), (-1, xj, "xt")], 0, zero)
    u_ = word("uu")
    L.mul_eq("u", z, xj, "z", "xt", [(1, u_, "uu")], 0)
    cc = word("cc")
    L.mul_eq("c", v, u_, "v", "uu", [(1, cc, "cc")], -2 * Bc)
    u2 = word("u2")
    L.mul_eq("u2e", u_, u_, "uu", "uu", [(1, u2, "u2")], 0)
    g = word("g")
    L.mul_eq("ge", dd, zn, "dd", "zn", [(1, g, "g")], 0)
    gg = word("gg")
    L.mul_eq("gge", g, g, "g", "g", [(1, gg, "gg")], 0)
    hh = word("hhe")
    L.mul_eq("hhe", cc, zn, "cc", "zn", [(1, hh, "hhe")], 0)
    poly = {}
    for t, vv in enumerate(gg):
        poly[(vv,)] = poly.get((vv,), 0) + (1 << t)
    for t, vv in enumerate(hh):
        poly[(vv,)] = poly.get((vv,), 0) - 2 * (1 << t)
    for t, vv in enumerate(u2):
        poly[(vv,)] = poly.get((vv,), 0) + (1 << t)
    for t, vv in enumerate(v):
        poly[(vv,)] = poly.get((vv,), 0) - 4 * Bc * (1 << t)
    L.congruent(poly, 0, "s3", zero)
    Q.finalize()
    st = Q.stats()
    st['vars'] -= s              # z belongs to the previous step
    return st


def modmul(s, mode='binary', kdepth=0, kmin=8, chunk=16, seed=1, wa=None, wb=None, toom=0):
    """cost of a single  A*B == W (mod p)  with |A| = wa, |B| = wb bits."""
    rnd = random.Random(seed)
    p = _prime(s, rnd)
    wa = wa or s
    wb = wb or s
    L = Ladder2(p, chunk=chunk, mode=mode, kdepth=kdepth, kmin=kmin, toom=toom)
    Q = L.qb
    zero = lambda wv: 0
    A = Q.word("A", wa, zero)
    Bw = Q.word("B", wb, zero)
    W = Q.word("W", s, zero)
    L.mul_eq("mm", A, Bw, "A", "B", [(1, W, "W")], 0)
    Q.finalize()
    st = Q.stats()
    st['vars'] -= (wa + wb)      # the operands are somebody else's
    return st


if __name__ == '__main__':
    import json
    import time
    s = 256
    print(f"{'what':>34} {'vars':>10} {'couplers':>13} {'AND':>9} {'|J|':>6}")
    rows = {}
    for kd in (0, 1, 2, 3, 4, 5, 6):
        t0 = time.time()
        st = modmul(s, kdepth=kd)
        rows[f'modmul_k{kd}'] = st
        print(f"{'modmul 256x256 kdepth=%d' % kd:>34} {st['vars']:10,d} "
              f"{st['couplers']:13,d} {st['and_vars']:9,d} "
              f"2^{st['dynamic_range_bits']:<4d} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({k: dict(v) for k, v in rows.items()},
              open('modmul256.json', 'w'), indent=1)
