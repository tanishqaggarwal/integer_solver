#!/usr/bin/env python3
"""ladder2.py -- the windowed comb ladder, rebuilt on the squeezed modmul.

Same encoding as ../ladder.py:build_win (one-hot comb windows, affine addition
with lam*d = e plus an explicit d != 0 gadget, sequential-counter one-hot), with
four changes, each measured separately in FINDINGS.md:

  * every modular reduction writes p in non-adjacent form, so a quotient bit
    contributes 6 column entries instead of popcount(p) = 250;
  * every table constant is written in NAF too, ~85 entries instead of ~128;
  * the three multiplications per group addition are Karatsuba / Toom trees;
  * lam*lam is a squaring (it already was -- the AND cache collapses it).

Faithfulness is re-established from scratch by demo_win2.py, which enumerates
every candidate scalar on a small curve exactly as demo_win.py does.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmqb import MMQB, best_split, bin_split                     # noqa: E402
from mm import build_modmul                                      # noqa: E402


def _qword(Q, p, name, lo, hi, valfn):
    """quotient word for an identity whose left side lies in [lo, hi]."""
    qlo, qhi = lo // p, hi // p
    nb = max(0, (qhi - qlo).bit_length())
    q = Q.mkword(f"q:{name}", nb, lambda wv, f=valfn, p=p, qlo=qlo: f(wv) // p - qlo)
    return q, qlo


def lin_word(Q, p, name, terms, const, fn, naf=True):
    """fresh s-bit word  ==  const + sum coef*word   (mod p)."""
    s = p.bit_length()
    W = Q.mkword(name, s, fn)
    split = best_split if naf else bin_split
    lo = hi = const
    for coef, V in terms:
        (hi, lo) = ((hi + coef * ((1 << len(V)) - 1), lo) if coef >= 0
                    else (hi, lo + coef * ((1 << len(V)) - 1)))
    lo -= (1 << s) - 1
    q, qlo = _qword(Q, p, name, lo, hi,
                    lambda wv, terms=terms, const=const, W=W:
                    const + sum(c * V.val(wv) for c, V in terms) - W.val(wv))
    tl = []
    for coef, V in terms:
        for sg, sh in split(coef):
            for t, b in enumerate(V.bits):
                tl.append(((b,), sg, sh + t))
    for t, b in enumerate(W.bits):
        tl.append(((b,), -1, t))
    for t, b in enumerate(q.bits):
        for sg, sh in split(p):
            tl.append(((b,), -sg, sh + t))
    Q.assert_terms(tl, bin_split(const) + bin_split(-p * qlo), f"lin:{name}")
    return W


def eq_const(Q, p, W, tgt, tag, naf=True):
    """assert  W == tgt  (mod p)."""
    split = best_split if naf else bin_split
    q, qlo = _qword(Q, p, tag, -tgt, (1 << len(W)) - 1 - tgt,
                    lambda wv, W=W, tgt=tgt: W.val(wv) - tgt)
    tl = [((b,), 1, t) for t, b in enumerate(W.bits)]
    for t, b in enumerate(q.bits):
        for sg, sh in split(p):
            tl.append(((b,), -sg, sh + t))
    Q.assert_terms(tl, bin_split(-tgt) + bin_split(-p * qlo), tag)


def not_equal(Q, p, W, c, tag):
    """assert the word W != the constant c, via hamming(W, c) = 1 + slack."""
    s = len(W)
    nb = max(1, (s - 1).bit_length())
    slack = Q.mkword(f"ne:{tag}", nb,
                     lambda wv, W=W, c=c, s=s:
                     bin((W.val(wv) ^ c) & ((1 << s) - 1)).count('1') - 1)
    tl, const = [], -1
    for i, v in enumerate(W.bits):
        if (c >> i) & 1:
            tl.append(((v,), -1, 0))
            const += 1
        else:
            tl.append(((v,), 1, 0))
    for t, v in enumerate(slack.bits):
        tl.append(((v,), -1, t))
    Q.assert_terms(tl, bin_split(const), tag)


def build_win2(p, table, T, w, mode='wallace', mult='karatsuba', leaf=32,
               red='naf', naf_consts=True, chunk=16, verbose=False):
    """table[j][t] = (t+1)*(2^{wj} G).  Constraint: sum_j table[j][d_j] == T."""
    M = len(table)
    D = 1 << w
    s = p.bit_length()
    Q = MMQB(chunk=chunk, mode=mode)
    split = best_split if naf_consts else bin_split

    U = []
    for j in range(M):
        u = []
        for t in range(D):
            v = Q.new(f"u{j}_{t}", 'input')
            Q.trace.append(('word', f"u{j}_{t}", [v],
                            (lambda wv, j=j, t=t: 1 if wv[f"_u{j}"] == t else 0)))
            u.append(v)
        prev = None
        for t, v in enumerate(u[:-1]):                 # sequential-counter one-hot
            pv = Q.new(f"p{j}_{t}", 'onehot')
            Q.trace.append(('word', f"p{j}_{t}", [pv],
                            (lambda wv, j=j, t=t: 1 if wv[f"_u{j}"] <= t else 0)))
            lin = {pv: -1, v: 1}
            if prev is not None:
                lin[prev] = 1
            Q.add_square(lin, 0)
            prev = pv
        Q.add_square({prev: 1, u[-1]: 1}, -1)
        U.append(u)

    def sel(name, j, coord):
        """word == sum_t u[j][t] * table[j][t][coord]  (mod p), constants in NAF."""
        W = Q.mkword(name, s,
                     lambda wv, j=j, coord=coord: table[j][wv[f"_u{j}"]][coord] % p)
        hi = max(table[j][t][coord] % p for t in range(D))
        q, qlo = _qword(Q, p, name, -((1 << s) - 1), hi,
                        lambda wv, j=j, coord=coord, W=W:
                        table[j][wv[f"_u{j}"]][coord] % p - W.val(wv))
        tl = []
        for t in range(D):
            for sg, sh in split(table[j][t][coord] % p):
                tl.append(((U[j][t],), sg, sh))
        for t, b in enumerate(W.bits):
            tl.append(((b,), -1, t))
        for t, b in enumerate(q.bits):
            for sg, sh in split(p):
                tl.append(((b,), -sg, sh + t))
        Q.assert_terms(tl, bin_split(-p * qlo), f"sel:{name}")
        return W

    x1 = sel("x1_0", 0, 0)
    y1 = sel("y1_0", 0, 1)

    for j in range(1, M):
        x2 = sel(f"x2_{j}", j, 0)
        y2 = sel(f"y2_{j}", j, 1)
        X1, Y1 = x1, y1
        d = lin_word(Q, p, f"d{j}", [(1, x2), (-1, X1)], 0,
                     lambda wv, x2=x2, X1=X1: (x2.val(wv) - X1.val(wv)) % p,
                     naf_consts)
        e = lin_word(Q, p, f"e{j}", [(1, y2), (-1, Y1)], 0,
                     lambda wv, y2=y2, Y1=Y1: (y2.val(wv) - Y1.val(wv)) % p,
                     naf_consts)
        lam = Q.mkword(f"lam{j}", s,
                       lambda wv, d=d, e=e, p=p:
                       e.val(wv) * pow(d.val(wv), -1, p) % p)
        build_modmul(Q, p, lam, d, [(e, 1)], mult=mult, leaf=leaf, red=red,
                     tag=f"lam{j}")
        for c in (0, p):
            not_equal(Q, p, d, c, f"ne{j}_{c}")
        x3 = Q.mkword(f"x3_{j}", s,
                      lambda wv, lam=lam, X1=X1, x2=x2, p=p:
                      (lam.val(wv) ** 2 - X1.val(wv) - x2.val(wv)) % p)
        build_modmul(Q, p, lam, lam, [(x3, 1), (X1, 1), (x2, 1)], mult=mult,
                     leaf=leaf, red=red, tag=f"x3e{j}")
        mm = lin_word(Q, p, f"m{j}", [(1, X1), (-1, x3)], 0,
                      lambda wv, X1=X1, x3=x3: (X1.val(wv) - x3.val(wv)) % p,
                      naf_consts)
        y3 = Q.mkword(f"y3_{j}", s,
                      lambda wv, lam=lam, mm=mm, Y1=Y1, p=p:
                      (lam.val(wv) * mm.val(wv) - Y1.val(wv)) % p)
        build_modmul(Q, p, lam, mm, [(y3, 1), (Y1, 1)], mult=mult, leaf=leaf,
                     red=red, tag=f"y3e{j}")
        x1, y1 = x3, y3
        if verbose:
            print(f"  window {j}/{M}: vars={Q.n:,}")

    Tx, Ty = T
    eq_const(Q, p, x1, Tx % p, 'finalX', naf_consts)
    eq_const(Q, p, y1, Ty % p, 'finalY', naf_consts)
    Q.finalize()
    return Q, U
