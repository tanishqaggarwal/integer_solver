#!/usr/bin/env python3
"""resources.py -- measured cost model of the encoding, extrapolated to the real instance."""
import random, json, sys
from qubo import QB
from ladder import Ladder


def cost(fn, mode):
    L = Ladder(P_DUMMY, mode=mode)
    fn(L)
    L.qb.finalize()
    st = L.qb.stats()
    return st['vars'], st['couplers']


def marginal_window(s, w, mode, seed=1, neq=True):
    """exact marginal cost of one comb window: 2 table look-ups + 4 modmuls + 4 linear words."""
    rnd = random.Random(seed)
    p = (1 << (s - 1)) + 2 * rnd.randrange(1 << (s - 3)) + 1
    D = 1 << w
    L = Ladder(p, mode=mode)
    Q = L.qb
    consts = [rnd.randrange(p) for _ in range(2 * D)]

    def word(name): return Q.word(name, s, lambda wv: 0)

    def lin(name, terms, const):
        wd = word(name)
        poly = {}
        for coef, bits in terms:
            for t, v in enumerate(bits): poly[(v,)] = poly.get((v,), 0) + coef * (1 << t)
        for t, v in enumerate(wd): poly[(v,)] = poly.get((v,), 0) - (1 << t)
        L.congruent(poly, const, f"lin:{name}", lambda wv: 0)
        return wd

    # --- the two table look-ups (one-hot MUX) ---
    sel = [Q.new(f"u{t}", 'input') for t in range(D)]
    Q.add_square({v: 1 for v in sel}, -1)
    x2 = lin("x2", [(consts[t], [sel[t]]) for t in range(D)], 0)
    y2 = lin("y2", [(consts[D + t], [sel[t]]) for t in range(D)], 0)
    # --- previous accumulator ---
    x1, y1 = word("x1"), word("y1")
    d = lin("d", [(1, x2), (-1, x1)], 0)
    e = lin("e", [(1, y2), (-1, y1)], 0)
    lam = word("lam")
    if neq:
        from ladder import not_equal
        L.mul_eq("lam", lam, d, "lam", "d", [(1, e, "e")], 0)
        for c in (0, p): not_equal(L, Q, d, "d", c, f"ne{c}")
    else:
        dinv = word("dinv")
        L.mul_eq("inv", d, dinv, "d", "dinv", [], 1)
        L.mul_eq("lam", e, dinv, "e", "dinv", [(1, lam, "lam")], 0)
    x3 = word("x3")
    L.mul_eq("x3", lam, lam, "lam", "lam", [(1, x3, "x3"), (1, x1, "x1"), (1, x2, "x2")], 0)
    mm = lin("m", [(1, x1), (-1, x3)], 0)
    y3 = word("y3")
    L.mul_eq("y3", lam, mm, "lam", "m", [(1, y3, "y3"), (1, y1, "y1")], 0)
    Q.finalize()
    st = Q.stats()
    # subtract the two carry-in registers x1,y1 that really belong to the previous window
    return st['vars'] - 2 * s, st['couplers'], st['dynamic_range_bits']


def fit_quadratic(xs, ys):
    """least squares y = a x^2 + b x + c"""
    import itertools
    n = len(xs)
    A = [[x*x, x, 1.0] for x in xs]
    AtA = [[sum(A[k][i]*A[k][j] for k in range(n)) for j in range(3)] for i in range(3)]
    Aty = [sum(A[k][i]*ys[k] for k in range(n)) for i in range(3)]
    for i in range(3):
        pv = max(range(i, 3), key=lambda r: abs(AtA[r][i]))
        AtA[i], AtA[pv] = AtA[pv], AtA[i]; Aty[i], Aty[pv] = Aty[pv], Aty[i]
        f = AtA[i][i]
        AtA[i] = [v/f for v in AtA[i]]; Aty[i] /= f
        for r in range(3):
            if r != i:
                g = AtA[r][i]
                AtA[r] = [AtA[r][c]-g*AtA[i][c] for c in range(3)]; Aty[r] -= g*Aty[i]
    return Aty


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'wallace'
    print(f"=== measured marginal cost of ONE comb window (mode={mode}) ===")
    print(f"{'s':>5} {'w':>3} {'vars':>12} {'couplers':>12} {'|J| bits':>9}")
    data = {}
    for w in (1, 4, 8):
        for s in (16, 24, 32, 48, 64):
            v, c, jr = marginal_window(s, w, mode)
            data[(s, w)] = (v, c)
            print(f"{s:5d} {w:3d} {v:12d} {c:12d} {jr:9d}")
    # model:  vars(s,w) = a s^2 + b s + c + mux(w) * s-ish   -> fit per w, then in w
    fits = {}
    for w in (1, 4, 8):
        xs = [16, 24, 32, 48, 64]
        fits[w] = fit_quadratic(xs, [data[(s, w)][0] for s in xs])
        a, b, c = fits[w]
        print(f"  fit w={w}: vars ~ {a:.2f} s^2 + {b:.1f} s + {c:.0f}")
    json.dump({f"{s},{w}": data[(s, w)] for (s, w) in data}, open(f'cost_{mode}.json', 'w'))
