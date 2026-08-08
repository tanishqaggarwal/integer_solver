#!/usr/bin/env python3
"""ladder.py -- build the QUBO for  sum_i b_i * (2^i G) == T  on y^2 = x^3 + B / F_p."""
from qubo import QB


class Ladder:
    def __init__(self, p, chunk=16, mode='binary'):
        self.p = p
        self.s = p.bit_length()
        self.qb = QB(chunk=chunk, mode=mode)

    # ---- assert  poly + const == 0  (mod p)  with an explicit quotient word ----
    def congruent(self, poly, const, tag, valfn):
        lo = const + sum(min(0, c) for c in poly.values())
        hi = const + sum(max(0, c) for c in poly.values())
        qlo, qhi = lo // self.p, hi // self.p
        nb = max(0, (qhi - qlo).bit_length())
        qbits = self.qb.word(f"q:{tag}", nb, lambda wv, f=valfn, b=qlo: f(wv) // self.p - b)
        poly = dict(poly)
        for t, v in enumerate(qbits):
            poly[(v,)] = poly.get((v,), 0) - self.p * (1 << t)
        self.qb.assert_zero(poly, const - self.p * qlo, tag)

    # ---- W == A*B (mod p) ----
    def mul_word(self, name, A, B, nameA, nameB):
        W = self.qb.word(name, self.s,
                         lambda wv, a=nameA, b=nameB: wv[a] * wv[b] % self.p)
        poly = {}
        for i, u in enumerate(A):
            for j, v in enumerate(B):
                m = (u, v) if u <= v else (v, u)
                poly[m] = poly.get(m, 0) + (1 << (i + j))
        for t, v in enumerate(W):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        self.congruent(poly, 0, f"mul:{name}",
                       lambda wv, a=nameA, b=nameB, n=name: wv[a] * wv[b] - wv[n])
        return W

    # ---- A*B == (sum terms + const) (mod p)   [no new result word] ----
    def mul_eq(self, tag, A, B, nameA, nameB, terms, const):
        poly = {}
        for i, u in enumerate(A):
            for j, v in enumerate(B):
                m = (u, v) if u <= v else (v, u)
                poly[m] = poly.get(m, 0) + (1 << (i + j))
        for coef, bits, _nm in terms:
            for t, v in enumerate(bits):
                poly[(v,)] = poly.get((v,), 0) - coef * (1 << t)

        def val(wv, a=nameA, b=nameB, terms=terms, const=const):
            return wv[a] * wv[b] - const - sum(c * wv[nm] for c, _bits, nm in terms)
        self.congruent(poly, -const, f"eq:{tag}", val)


def build(p, B, pts, T, chunk=16, mode='binary', verbose=False):
    """pts[i] = (X_i, Y_i) = 2^i G  (affine, canonical).  T = target point.
       Returns (L, bvars, ord_names)."""
    m = len(pts)
    L = Ladder(p, chunk=chunk, mode=mode)
    Q = L.qb
    s = L.s

    bv = [Q.new(f"b{i}", 'input') for i in range(m)]
    for i, v in enumerate(bv):
        Q.trace.append(('word', f"b{i}", [v], (lambda wv, i=i: wv[f"_b{i}"])))

    def W(name, terms, const, fn):
        """terms: list of (coef, bitlist, wordname)"""
        w = Q.word(name, s, fn)
        poly = {}
        for coef, bits, _ in terms:
            for t, v in enumerate(bits):
                poly[(v,)] = poly.get((v,), 0) + coef * (1 << t)
        for t, v in enumerate(w):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        L.congruent(poly, const, f"lin:{name}",
                    lambda wv, terms=terms, const=const, name=name:
                    const + sum(c * wv[nm] for c, _b, nm in terms) - wv[name])
        return w

    X0, Y0 = pts[0]
    x1 = W("x1_0", [], X0, lambda wv: X0 % p)
    y1 = W("y1_0", [(2 * Y0, [bv[0]], "b0")], -Y0, lambda wv: (2 * wv["_b0"] - 1) * Y0 % p)
    nx, ny = "x1_0", "y1_0"

    for i in range(1, m):
        Xi, Yi = pts[i]
        d = W(f"d{i}", [(-1, x1, nx)], Xi, lambda wv, Xi=Xi, nx=nx: (Xi - wv[nx]) % p)
        e = W(f"e{i}", [(2 * Yi, [bv[i]], f"b{i}"), (-1, y1, ny)], -Yi,
              lambda wv, Yi=Yi, ny=ny, i=i: ((2 * wv[f"_b{i}"] - 1) * Yi - wv[ny]) % p)
        dinv = Q.word(f"dinv{i}", s, lambda wv, i=i: pow(wv[f"d{i}"], -1, p))
        L.mul_eq(f"inv{i}", d, dinv, f"d{i}", f"dinv{i}", [], 1)
        lam = L.mul_word(f"lam{i}", e, dinv, f"e{i}", f"dinv{i}")
        x3 = Q.word(f"x3_{i}", s,
                    lambda wv, i=i, Xi=Xi, nx=nx: (wv[f"lam{i}"] ** 2 - wv[nx] - Xi) % p)
        L.mul_eq(f"x3e{i}", lam, lam, f"lam{i}", f"lam{i}",
                 [(1, x3, f"x3_{i}"), (1, x1, nx)], Xi)
        mm = W(f"m{i}", [(1, x1, nx), (-1, x3, f"x3_{i}")], 0,
               lambda wv, i=i, nx=nx: (wv[nx] - wv[f"x3_{i}"]) % p)
        y3 = Q.word(f"y3_{i}", s,
                    lambda wv, i=i, ny=ny: (wv[f"lam{i}"] * wv[f"m{i}"] - wv[ny]) % p)
        L.mul_eq(f"y3e{i}", lam, mm, f"lam{i}", f"m{i}",
                 [(1, y3, f"y3_{i}"), (1, y1, ny)], 0)
        x1, y1, nx, ny = x3, y3, f"x3_{i}", f"y3_{i}"
        if verbose and i % 32 == 0:
            print(f"  step {i}: vars={Q.n}")

    # final: running point == T  (mod p)
    Tx, Ty = T
    L.congruent({(v,): (1 << t) for t, v in enumerate(x1)}, -Tx, "finalX",
                lambda wv, nx=nx: wv[nx] - Tx)
    L.congruent({(v,): (1 << t) for t, v in enumerate(y1)}, -Ty, "finalY",
                lambda wv, ny=ny: wv[ny] - Ty)
    Q.finalize()
    return L, bv


def build_win(p, B, table, T, w, chunk=16, mode='binary', neq=False, verbose=False):
    """Windowed comb encoding.

    table[j][t] = (t+1) * (2^{w j} G)  for t = 0 .. 2^w-1   (window j, digit t+1)
    Constraint:  sum_j table[j][dig_j] == T     with exactly one digit per window.
    The caller has already folded the +1 digit offset into T.
    """
    M = len(table)
    D = 1 << w
    L = Ladder(p, chunk=chunk, mode=mode)
    Q = L.qb
    s = L.s

    # one-hot digit selectors
    U = []
    for j in range(M):
        u = [Q.new(f"u{j}_{t}", 'input') for t in range(D)]
        for t, v in enumerate(u):
            Q.trace.append(('word', f"u{j}_{t}", [v], (lambda wv, j=j, t=t: wv[f"_u{j}"] == t and 1 or 0)))
        # "exactly one digit", encoded as a sequential counter rather than
        # (sum u - 1)^2.  The naive square makes a K_D -- with D = 2^w that clique
        # alone dominates minor-embedding.  The prefix chain p_t = p_{t-1} + u_t
        # says the same thing with 3-term penalties, so the largest clique in the
        # whole Hamiltonian drops to the arithmetic's own width.
        prev = None
        for t, v in enumerate(u[:-1]):
            pv = Q.new(f"p{j}_{t}", 'onehot')
            Q.trace.append(('word', f"p{j}_{t}", [pv],
                            (lambda wv, j=j, t=t: 1 if wv[f"_u{j}"] <= t else 0)))
            lin = {pv: -1, v: 1}
            if prev is not None: lin[prev] = 1
            Q.add_square(lin, 0)
            prev = pv
        Q.add_square({prev: 1, u[-1]: 1}, -1)        # the last prefix must reach 1
        U.append(u)

    def W(name, terms, const, fn):
        wd = Q.word(name, s, fn)
        poly = {}
        for coef, bits, _ in terms:
            for t, v in enumerate(bits):
                poly[(v,)] = poly.get((v,), 0) + coef * (1 << t)
        for t, v in enumerate(wd):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        L.congruent(poly, const, f"lin:{name}",
                    lambda wv, terms=terms, const=const, name=name:
                    const + sum(c * wv[nm] for c, _b, nm in terms) - wv[name])
        return wd

    def select(name, j, coord):
        """word == sum_t u[j][t] * table[j][t][coord]   (mod p)"""
        terms = [(table[j][t][coord], [U[j][t]], f"u{j}_{t}") for t in range(D)]
        return W(name, terms, 0, lambda wv, j=j, coord=coord: table[j][wv[f"_u{j}"]][coord] % p)

    x1 = select("x1_0", 0, 0); nx = "x1_0"
    y1 = select("y1_0", 0, 1); ny = "y1_0"

    for j in range(1, M):
        x2 = select(f"x2_{j}", j, 0); n2x = f"x2_{j}"
        y2 = select(f"y2_{j}", j, 1); n2y = f"y2_{j}"
        d = W(f"d{j}", [(1, x2, n2x), (-1, x1, nx)], 0,
              lambda wv, n2x=n2x, nx=nx: (wv[n2x] - wv[nx]) % p)
        e = W(f"e{j}", [(1, y2, n2y), (-1, y1, ny)], 0,
              lambda wv, n2y=n2y, ny=ny: (wv[n2y] - wv[ny]) % p)
        if neq:
            # 3-multiplication form: lam*d == e, with d != 0 enforced directly.
            # d is a word in [0,2^s); d == 0 (mod p) iff d in {0, p}, so two
            # "not equal to this constant" gadgets close the degenerate-division
            # loophole for ~2*log2(s) ancillas instead of a whole modular inverse.
            lam = Q.word(f"lam{j}", s, lambda wv, j=j:
                         wv[f"e{j}"] * pow(wv[f"d{j}"], -1, p) % p)
            L.mul_eq(f"lam{j}", lam, d, f"lam{j}", f"d{j}", [(1, e, f"e{j}")], 0)
            for c in (0, p):
                not_equal(L, Q, d, f"d{j}", c, f"ne{j}_{c}")
        else:
            dinv = Q.word(f"dinv{j}", s, lambda wv, j=j: pow(wv[f"d{j}"], -1, p))
            L.mul_eq(f"inv{j}", d, dinv, f"d{j}", f"dinv{j}", [], 1)
            lam = L.mul_word(f"lam{j}", e, dinv, f"e{j}", f"dinv{j}")
        x3 = Q.word(f"x3_{j}", s, lambda wv, j=j, nx=nx, n2x=n2x:
                    (wv[f"lam{j}"] ** 2 - wv[nx] - wv[n2x]) % p)
        L.mul_eq(f"x3e{j}", lam, lam, f"lam{j}", f"lam{j}",
                 [(1, x3, f"x3_{j}"), (1, x1, nx), (1, x2, n2x)], 0)
        mm = W(f"m{j}", [(1, x1, nx), (-1, x3, f"x3_{j}")], 0,
               lambda wv, j=j, nx=nx: (wv[nx] - wv[f"x3_{j}"]) % p)
        y3 = Q.word(f"y3_{j}", s, lambda wv, j=j, ny=ny:
                    (wv[f"lam{j}"] * wv[f"m{j}"] - wv[ny]) % p)
        L.mul_eq(f"y3e{j}", lam, mm, f"lam{j}", f"m{j}",
                 [(1, y3, f"y3_{j}"), (1, y1, ny)], 0)
        x1, y1, nx, ny = x3, y3, f"x3_{j}", f"y3_{j}"
        if verbose: print(f"  window {j}/{M}: vars={Q.n}")

    Tx, Ty = T
    L.congruent({(v,): (1 << t) for t, v in enumerate(x1)}, -Tx, "finalX",
                lambda wv, nx=nx: wv[nx] - Tx)
    L.congruent({(v,): (1 << t) for t, v in enumerate(y1)}, -Ty, "finalY",
                lambda wv, ny=ny: wv[ny] - Ty)
    Q.finalize()
    return L, U


def not_equal(L, Q, word, wname, c, tag):
    """assert the s-bit word != the constant c, via  hamming(word, c) = 1 + slack."""
    s = len(word)
    nb = max(1, (s - 1).bit_length())
    slack = Q.word(f"ne:{tag}", nb, lambda wv, wname=wname, c=c, s=s:
                   bin((wv[wname] ^ c) & ((1 << s) - 1)).count('1') - 1)
    poly, const = {}, -1
    for i, v in enumerate(word):
        if (c >> i) & 1: poly[(v,)] = poly.get((v,), 0) - 1; const += 1
        else: poly[(v,)] = poly.get((v,), 0) + 1
    for t, v in enumerate(slack): poly[(v,)] = poly.get((v,), 0) - (1 << t)
    Q.assert_zero(poly, const, tag)
