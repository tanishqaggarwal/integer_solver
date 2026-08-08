#!/usr/bin/env python3
"""mm.py -- ONE modular multiplication, every way we know how to build it.

    build_modmul(Q, p, A, B, C, ...)     asserts   A*B == C  (mod p)

Options, all measurable independently:

  mult    'schoolbook' | 'karatsuba' | 'toom3'      partial-product structure
  leaf    recursion leaf width (karatsuba / toom3)
  red     'quotient'   general quotient word, coefficient p binary-expanded
          'naf'        same quotient word, coefficient p in non-adjacent form
                       (for p = 2^256-2^32-977 that is 6 terms, not 250)
          'fold'       explicit pseudo-Mersenne folding  X_lo + c*X_hi
  square  A is B: the AND cache collapses a_i b_j / a_j b_i, s(s+1)/2 not s^2

The product is never materialised: a multiplier returns a LIST OF SIGNED
POWER-OF-TWO TERMS which is spliced straight into the reduction identity.
Recursive multipliers do materialise their sub-products, because the
recombination  A*B = P0 + M*2^h + P1*2^2h  needs them as words.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmqb import MMQB, Wrd, best_split, bin_split                # noqa: E402


# ------------------------------------------------------------- multipliers
def school_terms(A, B, shift=0):
    """terms for A*B, s_A * s_B AND monomials."""
    out = []
    for i, u in enumerate(A.bits):
        for j, v in enumerate(B.bits):
            out.append(((u, v), 1, i + j + shift))
    return out


def _add_word(Q, name, parts, nbits, tag):
    """fresh word == sum_k sign_k * 2^shift_k * word_k ; asserted over Z."""
    def fn(wv, parts=parts):
        return sum(sg * (v.val(wv) << sh) for v, sg, sh in parts)
    w = Q.mkword(name, nbits, fn)
    terms = []
    for v, sg, sh in parts:
        for t, b in enumerate(v.bits):
            terms.append(((b,), sg, t + sh))
    for t, b in enumerate(w.bits):
        terms.append(((b,), -1, t))
    Q.assert_terms(terms, [], tag)
    return w


def _mul_word(Q, name, A, B, mult, leaf, tag, counter):
    """fresh word == A*B (exact, over Z)."""
    nb = len(A) + len(B)
    w = Q.mkword(name, nb, lambda wv, A=A, B=B: A.val(wv) * B.val(wv))
    terms = product_terms(Q, A, B, mult, leaf, tag, counter)
    for t, b in enumerate(w.bits):
        terms.append(((b,), -1, t))
    Q.assert_terms(terms, [], tag)
    return w


def product_terms(Q, A, B, mult='schoolbook', leaf=32, tag='m', counter=None):
    """signed power-of-two terms whose value is exactly A*B."""
    if counter is None:
        counter = [0]
    n = max(len(A), len(B))
    if mult == 'schoolbook' or n <= leaf or min(len(A), len(B)) <= 2:
        return school_terms(A, B)
    counter[0] += 1
    tg = f"{tag}.{counter[0]}"
    if mult == 'karatsuba':
        return _kara(Q, A, B, leaf, tg, counter)
    if mult == 'toom3':
        # Toom-3's evaluation words are h+3 bits wide; recursion only pays if
        # that is strictly narrower than the input, otherwise fall back.
        h = (n + 2) // 3
        if h + 4 >= n:
            counter[0] -= 1
            return school_terms(A, B)
        return _toom3(Q, A, B, leaf, tg, counter)
    raise ValueError(mult)


# ----------------------------------------------------------- Karatsuba
def _kara(Q, A, B, leaf, tag, counter):
    n = max(len(A), len(B))
    h = n // 2
    A0, A1 = A.slice(0, min(h, len(A))), A.slice(h, len(A)) if len(A) > h else _zero()
    B0, B1 = B.slice(0, min(h, len(B))), B.slice(h, len(B)) if len(B) > h else _zero()
    same = A.bits == B.bits
    U = _add_word(Q, f"ku{tag}", [(A0, 1, 0), (A1, 1, 0)],
                  max(len(A0), len(A1)) + 1, f"ku{tag}")
    V = U if same else _add_word(Q, f"kv{tag}", [(B0, 1, 0), (B1, 1, 0)],
                                 max(len(B0), len(B1)) + 1, f"kv{tag}")
    P0 = _mul_word(Q, f"kp0{tag}", A0, B0, 'karatsuba', leaf, f"kp0{tag}", counter)
    P1 = _mul_word(Q, f"kp1{tag}", A1, B1, 'karatsuba', leaf, f"kp1{tag}", counter)
    P2 = _mul_word(Q, f"kp2{tag}", U, V, 'karatsuba', leaf, f"kp2{tag}", counter)
    # Mid = P2 - P0 - P1 >= 0, materialised so every coefficient stays a single
    # power of two (folding the subtraction into a coefficient would turn one
    # term into h column entries).
    Mid = Q.mkword(f"kmid{tag}", len(P2),
                   lambda wv, P0=P0, P1=P1, P2=P2: P2.val(wv) - P0.val(wv) - P1.val(wv))
    terms = []
    for w, sg in ((P0, 1), (P1, 1), (Mid, 1), (P2, -1)):
        for t, b in enumerate(w.bits):
            terms.append(((b,), sg, t))
    Q.assert_terms(terms, [], f"kmid{tag}")
    out = []
    for w, sh in ((P0, 0), (Mid, h), (P1, 2 * h)):
        for t, b in enumerate(w.bits):
            out.append(((b,), 1, t + sh))
    return out


def _zero():
    return Wrd([], lambda wv: 0)


# ----------------------------------------------------------- Toom-Cook 3
# evaluation points 0, 1, 2, 3, infinity -- all non-negative, so every word in
# the tree stays unsigned and the existing column balancer is reused verbatim.
def _toom3(Q, A, B, leaf, tag, counter):
    n = max(len(A), len(B))
    h = (n + 2) // 3
    def part(W, k):
        lo, hi = k * h, min((k + 1) * h, len(W))
        return W.slice(lo, hi) if lo < len(W) else _zero()
    a = [part(A, k) for k in range(3)]
    b = [part(B, k) for k in range(3)]
    same = A.bits == B.bits
    # a(x) = a0 + a1 x + a2 x^2  at x = 1,2,3
    def ev(ws, pt, nm):
        parts = [(ws[k], 1, 0) for k in range(3) if len(ws[k])]
        coeffs = [pt ** k for k in range(3)]
        nbits = max(len(ws[k]) + coeffs[k].bit_length() for k in range(3) if len(ws[k])) + 2
        prts = []
        for k in range(3):
            if not len(ws[k]):
                continue
            for sg, sh in bin_split(coeffs[k]):
                prts.append((ws[k], sg, sh))
        return _add_word(Q, nm, prts, nbits, nm)
    av = {pt: ev(a, pt, f"ta{tag}_{pt}") for pt in (1, 2, 3)}
    bv = av if same else {pt: ev(b, pt, f"tb{tag}_{pt}") for pt in (1, 2, 3)}
    W = {}
    W[0] = _mul_word(Q, f"tw0{tag}", a[0], b[0], 'toom3', leaf, f"tw0{tag}", counter)
    W[4] = _mul_word(Q, f"tw4{tag}", a[2], b[2], 'toom3', leaf, f"tw4{tag}", counter)
    for pt in (1, 2, 3):
        W[pt + 10] = _mul_word(Q, f"tw{pt}{tag}", av[pt], bv[pt], 'toom3', leaf,
                               f"tw{pt}{tag}", counter)
    # C(x) = c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4 ; c0 = W0, c4 = W4.
    # W(1) = c0+c1+c2+c3+c4 ; W(2) = c0+2c1+4c2+8c3+16c4 ; W(3) = c0+3c1+9c2+27c3+81c4
    def cval(k):
        def f(wv, k=k):
            w0, w4 = W[0].val(wv), W[4].val(wv)
            y1 = W[11].val(wv) - w0 - w4
            y2 = W[12].val(wv) - w0 - 16 * w4
            y3 = W[13].val(wv) - w0 - 81 * w4
            # solve [1 1 1; 2 4 8; 3 9 27] [c1 c2 c3]^T = [y1 y2 y3]^T
            c3 = (y3 - 3 * y2 + 3 * y1) // 6
            c2 = (y2 - 2 * y1 - 6 * c3) // 2
            c1 = y1 - c2 - c3
            return (c1, c2, c3)[k - 1]
        return f
    nb = len(A) + len(B)
    c1 = Q.mkword(f"tc1{tag}", nb, cval(1))
    c2 = Q.mkword(f"tc2{tag}", nb, cval(2))
    c3 = Q.mkword(f"tc3{tag}", nb, cval(3))
    for pt, wk in ((1, 11), (2, 12), (3, 13)):
        terms = []
        for t, bb in enumerate(W[wk].bits):
            terms.append(((bb,), -1, t))
        for k, cw in ((0, W[0]), (1, c1), (2, c2), (3, c3), (4, W[4])):
            for sg, sh in bin_split(pt ** k):
                for t, bb in enumerate(cw.bits):
                    terms.append(((bb,), sg, t + sh))
        Q.assert_terms(terms, [], f"ti{tag}_{pt}")
    out = []
    for k, cw in ((0, W[0]), (1, c1), (2, c2), (3, c3), (4, W[4])):
        for t, bb in enumerate(cw.bits):
            out.append(((bb,), 1, t + k * h))
    return out


# ------------------------------------------------------------- reduction
def build_modmul(Q, p, A, B, C, mult='schoolbook', leaf=32, red='naf', tag='mm'):
    """assert A*B == C (mod p).  A, B, C are Wrd.  Returns nothing."""
    counter = [0]
    terms = product_terms(Q, A, B, mult, leaf, tag, counter)
    if red == 'fold':
        return _fold_reduce(Q, p, terms, len(A) + len(B), C, tag)
    split = best_split(p) if red == 'naf' else bin_split(p)
    hi = ((1 << len(A)) - 1) * ((1 << len(B)) - 1)
    lo = -((1 << len(C)) - 1)
    qlo, qhi = lo // p, hi // p
    nb = max(0, (qhi - qlo).bit_length())
    q = Q.mkword(f"q:{tag}", nb,
                 lambda wv, A=A, B=B, C=C, p=p, qlo=qlo:
                 (A.val(wv) * B.val(wv) - C.val(wv)) // p - qlo)
    for t, b in enumerate(C.bits):
        terms.append(((b,), -1, t))
    for t, b in enumerate(q.bits):
        for sg, sh in split:
            terms.append(((b,), -sg, sh + t))
    Q.assert_terms(terms, bin_split(-p * qlo), tag)


def _fold_reduce(Q, p, terms, nb_prod, C, tag):
    """explicit pseudo-Mersenne folding.  p = 2^k - c ; X = Xlo + 2^k Xhi
       ==> X == Xlo + c*Xhi (mod p).  Iterated until the residue fits, then a
       single small quotient closes it."""
    k = p.bit_length()
    c = (1 << k) - p
    assert 0 < c, "not pseudo-Mersenne"
    csplit = best_split(c)
    cur, curbits, stage = terms, nb_prod, 0
    while curbits > k + 2:
        stage += 1
        hi_bits = curbits - k
        lo = Q.mkword(f"fl{tag}{stage}", k, lambda wv: 0)   # values patched below
        hiw = Q.mkword(f"fh{tag}{stage}", hi_bits, lambda wv: 0)
        # value functions need the current polynomial's value; carry it along
        _bind_fold(Q, cur, lo, hiw, k)
        eq = list(cur)
        for t, b in enumerate(lo.bits):
            eq.append(((b,), -1, t))
        for t, b in enumerate(hiw.bits):
            eq.append(((b,), -1, t + k))
        Q.assert_terms(eq, [], f"{tag}:split{stage}")
        cur = [((b,), 1, t) for t, b in enumerate(lo.bits)]
        for t, b in enumerate(hiw.bits):
            for sg, sh in csplit:
                cur.append(((b,), sg, sh + t))
        curbits = max(k, hi_bits + c.bit_length()) + 1
    # final: cur - C = p*q with q tiny
    hi = (1 << curbits) - 1
    lo_ = -((1 << len(C)) - 1)
    qlo, qhi = lo_ // p, hi // p
    nbq = max(0, (qhi - qlo).bit_length())
    val_cur = _fold_value(Q, cur)
    q = Q.mkword(f"q:{tag}", nbq,
                 lambda wv, C=C, p=p, qlo=qlo, f=val_cur: (f(wv) - C.val(wv)) // p - qlo)
    eq = list(cur)
    for t, b in enumerate(C.bits):
        eq.append(((b,), -1, t))
    for t, b in enumerate(q.bits):
        for sg, sh in best_split(p):
            eq.append(((b,), -sg, sh + t))
    Q.assert_terms(eq, bin_split(-p * qlo), tag)


# The fold needs the *value* of a term list.  Terms are (monomial, sign, shift)
# over variables whose witness values are only known at replay time, so we build
# the evaluator out of the same trace machinery: every monomial here is a single
# bit of an already-named word, so we record (wordname, bitindex).
_BITOWNER = {}


def _register_bits(Q):
    """map variable id -> (word name, bit index) for every word built so far."""
    for op in Q.trace:
        if op[0] == 'word':
            _, name, bits, _fn = op
            for t, v in enumerate(bits):
                _BITOWNER[(id(Q), v)] = (name, t)


def _fold_value(Q, terms):
    _register_bits(Q)
    items = []
    for mono, sg, sh in terms:
        assert len(mono) == 1, "fold value needs linearised terms"
        items.append((_BITOWNER[(id(Q), mono[0])], sg, sh))

    def f(wv, items=items):
        tot = 0
        for (nm, t), sg, sh in items:
            tot += sg * (((wv[nm] >> t) & 1) << sh)
        return tot
    return f


def _bind_fold(Q, terms, lo, hiw, k):
    """patch the witness functions of lo/hi so they split the value of `terms`."""
    _register_bits(Q)
    items = []
    for mono, sg, sh in terms:
        if len(mono) == 1:
            items.append((('bit', _BITOWNER[(id(Q), mono[0])]), sg, sh))
        else:
            items.append((('and', tuple(_BITOWNER[(id(Q), v)] for v in mono)), sg, sh))

    def value(wv, items=items):
        tot = 0
        for kind, sg, sh in items:
            if kind[0] == 'bit':
                nm, t = kind[1]
                bit = (wv[nm] >> t) & 1
            else:
                bit = 1
                for nm, t in kind[1]:
                    bit &= (wv[nm] >> t) & 1
            tot += sg * (bit << sh)
        return tot
    for i, op in enumerate(Q.trace):
        if op[0] == 'word' and op[1] == lo.name:
            Q.trace[i] = ('word', lo.name, lo.bits,
                          lambda wv: value(wv) & ((1 << k) - 1))
        elif op[0] == 'word' and op[1] == hiw.name:
            Q.trace[i] = ('word', hiw.name, hiw.bits,
                          lambda wv: value(wv) >> k)
