#!/usr/bin/env python3
"""pieces.py -- measure the cost of every ATOMIC fragment of the encoding at s = 256.

Nothing here is extrapolated: each gadget is compiled for real with the same
`qubo.py` / `ladder.py` used by report.py, and the qubits belonging to its
*inputs* (words the previous stage already paid for) are subtracted, so what is
reported is the marginal cost of the fragment itself.

Fragments, from largest to smallest:

    win(w)   one comb window   = 2 selects + 1 EC addition           (= report.py's unit)
    add      one EC addition   = 2 lin + 3 modmul + 2 neq + 2 lin + 2 words
    mul      one 256x256 modular multiplication  a*b == c (mod p)
    sel(w)   one one-hot table look-up: word == sum_t u_t * const_t (mod p)
    lin      one linear congruence word:  word == (linear combo) (mod p)
    neq      one "word != constant" gadget
    col      one balanced column inside a modmul   (the true atom of the compiler)

Writes multirun/pieces256.json.
"""
import json, os, random, sys, time

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

from qubo import QB
from ladder import Ladder, not_equal

S = 256
# a 256-bit prime of the same shape as the instance's (p = 2^256 - 2^32 - 977)
P = (1 << 256) - (1 << 32) - 977


def _fresh(mode):
    L = Ladder(P, mode=mode)
    return L, L.qb


def _lin(L, name, terms, const=0):
    """word == const + sum coef*word_i  (mod p);  terms = [(coef, bitlist)]"""
    Q = L.qb
    wd = Q.word(name, S, lambda wv: 0)
    poly = {}
    for coef, bits in terms:
        for t, v in enumerate(bits):
            poly[(v,)] = poly.get((v,), 0) + coef * (1 << t)
    for t, v in enumerate(wd):
        poly[(v,)] = poly.get((v,), 0) - (1 << t)
    L.congruent(poly, const, f"lin:{name}", lambda wv: 0)
    return wd


def _stats(Q, sub_vars):
    Q.finalize()
    st = Q.stats()
    return dict(vars=st['vars'] - sub_vars, couplers=st['couplers'],
                jbits=st['dynamic_range_bits'])


# ---------------------------------------------------------------- fragments
def piece_mul(mode):
    """one 256x256 modular multiplication:  a*b == c (mod p), a,b,c pre-existing."""
    L, Q = _fresh(mode)
    a = Q.word("a", S, lambda wv: 0)
    b = Q.word("b", S, lambda wv: 0)
    c = Q.word("c", S, lambda wv: 0)
    L.mul_eq("mm", a, b, "a", "b", [(1, c, "c")], 0)
    return _stats(Q, 3 * S)


def piece_lin(mode):
    """one linear congruence word (2 operands)."""
    L, Q = _fresh(mode)
    a = Q.word("a", S, lambda wv: 0)
    b = Q.word("b", S, lambda wv: 0)
    _lin(L, "r", [(1, a), (-1, b)])
    return _stats(Q, 2 * S)


def piece_neq(mode):
    """one 'word != constant' gadget."""
    L, Q = _fresh(mode)
    a = Q.word("a", S, lambda wv: 0)
    not_equal(L, Q, a, "a", 0, "ne0")
    return _stats(Q, S)


def piece_sel(mode, w):
    """one one-hot table look-up producing one s-bit word (2^w entries)."""
    rnd = random.Random(7)
    L, Q = _fresh(mode)
    D = 1 << w
    u = [Q.new(f"u{t}", 'input') for t in range(D)]
    Q.add_square({v: 1 for v in u}, -1)
    consts = [rnd.randrange(P) for _ in range(D)]
    _lin(L, "x2", [(consts[t], [u[t]]) for t in range(D)])
    return _stats(Q, 0)


def piece_add(mode, neq=True):
    """one full EC addition (x1,y1,x2,y2 already exist)."""
    L, Q = _fresh(mode)
    x1 = Q.word("x1", S, lambda wv: 0)
    y1 = Q.word("y1", S, lambda wv: 0)
    x2 = Q.word("x2", S, lambda wv: 0)
    y2 = Q.word("y2", S, lambda wv: 0)
    d = _lin(L, "d", [(1, x2), (-1, x1)])
    e = _lin(L, "e", [(1, y2), (-1, y1)])
    lam = Q.word("lam", S, lambda wv: 0)
    if neq:
        L.mul_eq("lam", lam, d, "lam", "d", [(1, e, "e")], 0)
        for c in (0, P):
            not_equal(L, Q, d, "d", c, f"ne{c}")
    else:
        dinv = Q.word("dinv", S, lambda wv: 0)
        L.mul_eq("inv", d, dinv, "d", "dinv", [], 1)
        L.mul_eq("lam", e, dinv, "e", "dinv", [(1, lam, "lam")], 0)
    x3 = Q.word("x3", S, lambda wv: 0)
    L.mul_eq("x3", lam, lam, "lam", "lam", [(1, x3, "x3"), (1, x1, "x1"), (1, x2, "x2")], 0)
    m = _lin(L, "m", [(1, x1), (-1, x3)])
    y3 = Q.word("y3", S, lambda wv: 0)
    L.mul_eq("y3", lam, m, "lam", "m", [(1, y3, "y3"), (1, y1, "y1")], 0)
    return _stats(Q, 4 * S)


def piece_final(mode):
    """the terminal check: running point == T  (both coordinates)."""
    L, Q = _fresh(mode)
    x1 = Q.word("x1", S, lambda wv: 0)
    y1 = Q.word("y1", S, lambda wv: 0)
    L.congruent({(v,): (1 << t) for t, v in enumerate(x1)}, -12345, "finalX", lambda wv: 0)
    L.congruent({(v,): (1 << t) for t, v in enumerate(y1)}, -54321, "finalY", lambda wv: 0)
    return _stats(Q, 2 * S)


def piece_col(mode):
    """the compiler's true atom: ONE balanced column of a 256x256 modular multiply.

    Reported as (vars of a whole modmul) / (number of columns it balances).
    A modmul at s bits balances ~2s+1 columns.
    """
    mm = piece_mul(mode)
    ncol = 2 * S + 2                     # columns of the 2s-bit product plus quotient span
    return dict(vars=mm['vars'] / ncol, couplers=mm['couplers'] / ncol,
                jbits=mm['jbits'], ncol=ncol)


# ---------------------------------------------------------------- driver
if __name__ == '__main__':
    out = {}
    for mode in ('binary', 'wallace'):
        t0 = time.time()
        r = {}
        r['mul'] = piece_mul(mode)
        r['lin'] = piece_lin(mode)
        r['neq'] = piece_neq(mode)
        r['add'] = piece_add(mode, neq=True)
        r['add_inv'] = piece_add(mode, neq=False)
        r['final'] = piece_final(mode)
        r['col'] = piece_col(mode)
        for w in (1, 2, 4, 6, 8, 10, 12):
            r[f'sel{w}'] = piece_sel(mode, w)
        r['_secs'] = round(time.time() - t0, 1)
        out[mode] = r
        print(f"=== {mode}  ({r['_secs']}s) ===")
        print(f"{'fragment':>10} {'qubits':>12} {'couplers':>14} {'|J| bits':>9}")
        for k in ('col', 'neq', 'lin', 'final', 'sel1', 'sel2', 'sel4', 'sel6',
                  'sel8', 'sel10', 'sel12', 'mul', 'add_inv', 'add'):
            v = r[k]
            print(f"{k:>10} {v['vars']:12,.0f} {v['couplers']:14,.0f} {v['jbits']:9d}")
        print()
    json.dump(out, open(os.path.join(_HERE, 'pieces256.json'), 'w'), indent=1)
    print("wrote multirun/pieces256.json")
