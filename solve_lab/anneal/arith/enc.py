#!/usr/bin/env python3
"""enc.py -- alternative arithmetizations of

    find b_0..b_{m-1} in {0,1}  with  sum_i b_i (2^i G) == T   on  y^2 = x^3 + B / F_p

built on top of the existing qubo.QB / ladder.Ladder machinery (imported, never
modified).  Everything here is either

  * an EXACT re-encoding (ground states == solutions), verified exhaustively by
    arith/faith.py, or
  * a RELAXATION, in which case faith.py measures how many spurious ground
    states it actually admits.

Orthogonal knobs, measurable independently and in combination:

  mux=True      table look-ups compiled as a one-hot selector penalty instead of
                a weighted-sum congruence.  ZERO ancillas per table entry.
  kdepth=d      every s x s multiplication done by d levels of Karatsuba.
  signed=True   signed-digit comb: digits +-{1,3,...,2^w-1}, so the table holds
                2^{w-1} magnitudes and one sign bit flips y.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ladder import Ladder, not_equal      # noqa: E402
from qubo import QB                       # noqa: E402


# ===================================================================== Ladder2
class Ladder2(Ladder):
    """Ladder + Karatsuba multiplication + a zero-ancilla one-hot MUX."""

    def __init__(self, p, chunk=16, mode='binary', kdepth=0, kmin=8, toom=0):
        super().__init__(p, chunk=chunk, mode=mode)
        self.kdepth = kdepth
        self.kmin = kmin          # stop recursing at or below this operand width
        self.toom = toom          # how many TOP levels use Toom-3 instead of Karatsuba
        self._kc = 0

    # ---------------------------------------------------------- Karatsuba ----
    def _tag(self, base):
        self._kc += 1
        return f"{base}#{self._kc}"

    # ------------------------------------------------------------- Toom-3 ----
    def _toom3_poly(self, A, fA, B, fB, tag, depth, toom):
        """Toom-Cook 3-way: 5 sub-multiplies of width ~n/3 instead of 9.

        Evaluation points 0, 1, -1, 2, infinity.  A(-1) can be negative, so it is
        carried biased by 2^h and un-biased linearly.  The interpolation is NOT
        done by the usual divide-by-2-and-3 sequence: the base-2^h columns
        c_0..c_4 of the product are all non-negative, so they are materialised as
        words and the three Vandermonde rows are asserted directly as exact
        integer identities.  Same result, no exact-division gadgets."""
        Q = self.qb
        n = max(len(A), len(B))
        h = (n + 2) // 3
        sl = lambda X, i: X[i * h:(i + 1) * h] if i < 2 else X[2 * h:]
        pc = lambda f, i, hh=h: ((f(wv) >> (i * hh)) & ((1 << hh) - 1))

        def part(f, i):
            if i < 2:
                return lambda wv, f=f, i=i, hh=h: (f(wv) >> (i * hh)) & ((1 << hh) - 1)
            return lambda wv, f=f, hh=h: f(wv) >> (2 * hh)

        A0, A1, A2 = sl(A, 0), sl(A, 1), sl(A, 2)
        B0, B1, B2 = sl(B, 0), sl(B, 1), sl(B, 2)
        fa = [part(fA, i) for i in range(3)]
        fb = [part(fB, i) for i in range(3)]
        same = (A is B) or (A == B)

        def evalword(P, f, coefs, bias, nm):
            """word == bias + sum coefs[i]*part_i  (all non-negative)."""
            nbits = max(1, (bias + sum(c * ((1 << len(P[i])) - 1)
                                       for i, c in enumerate(coefs)
                                       if c > 0)).bit_length())
            W = Q.word(nm, nbits, lambda wv, f=f, c=coefs, b=bias:
                       b + sum(c[i] * f[i](wv) for i in range(3)))
            poly = {}
            for i, c in enumerate(coefs):
                for t, v in enumerate(P[i]):
                    poly[(v,)] = poly.get((v,), 0) + c * (1 << t)
            for t, v in enumerate(W):
                poly[(v,)] = poly.get((v,), 0) - (1 << t)
            Q.assert_zero(poly, bias, nm)
            return W, (lambda wv, nm=nm: wv[nm])

        PA, PB = [A0, A1, A2], [B0, B1, B2]
        Ap1, fAp1 = evalword(PA, fa, [1, 1, 1], 0, f"t1a:{tag}")
        Am1, fAm1 = evalword(PA, fa, [1, -1, 1], 1 << h, f"tma:{tag}")
        Ap2, fAp2 = evalword(PA, fa, [1, 2, 4], 0, f"t2a:{tag}")
        if same:
            Bp1, fBp1, Bm1, fBm1, Bp2, fBp2 = Ap1, fAp1, Am1, fAm1, Ap2, fAp2
        else:
            Bp1, fBp1 = evalword(PB, fb, [1, 1, 1], 0, f"t1b:{tag}")
            Bm1, fBm1 = evalword(PB, fb, [1, -1, 1], 1 << h, f"tmb:{tag}")
            Bp2, fBp2 = evalword(PB, fb, [1, 2, 4], 0, f"t2b:{tag}")

        P0, f0 = self._prod_word(A0, fa[0], B0, fb[0], tag + '0', depth, toom)
        P1, f1 = self._prod_word(Ap1, fAp1, Bp1, fBp1, tag + '1', depth, toom)
        PM, fM = self._prod_word(Am1, fAm1, Bm1, fBm1, tag + 'm', depth, toom)
        P2, f2 = self._prod_word(Ap2, fAp2, Bp2, fBp2, tag + '2', depth, toom)
        PI, fI = self._prod_word(A2, fa[2], B2, fb[2], tag + 'i', depth, toom)

        # pm1 = A(-1)B(-1) = PM - 2^h (Am1 + Bm1) + 2^{2h}
        pm1 = {}
        for t, v in enumerate(PM):
            pm1[(v,)] = pm1.get((v,), 0) + (1 << t)
        for W in ((Am1, Bm1) if not same else (Am1, Am1)):
            for t, v in enumerate(W):
                pm1[(v,)] = pm1.get((v,), 0) - (1 << (h + t))
        pm1_c = 1 << (2 * h)

        # c_i = base-2^h columns of the product; all >= 0
        def col(i):
            return lambda wv, i=i: sum(fa[j](wv) * fb[i - j](wv)
                                       for j in range(3) if 0 <= i - j < 3)
        nb = 2 * h + 3
        C = [None] * 5
        C[0], C[4] = P0, PI
        for i in (1, 2, 3):
            C[i] = Q.word(f"tc{i}:{tag}", nb, col(i))

        def lin(coefs):
            poly = {}
            for i, c in enumerate(coefs):
                for t, v in enumerate(C[i]):
                    poly[(v,)] = poly.get((v,), 0) + c * (1 << t)
            return poly

        def combine(a, b, const=0):
            poly = dict(a)
            for m, c in b.items():
                poly[m] = poly.get(m, 0) - c
            Q.assert_zero(poly, const, f"toomv:{tag}")

        combine(lin([1, 1, 1, 1, 1]), {(v,): (1 << t) for t, v in enumerate(P1)})
        combine(lin([1, -1, 1, -1, 1]), pm1, -pm1_c)
        combine(lin([1, 2, 4, 8, 16]), {(v,): (1 << t) for t, v in enumerate(P2)})

        out = {}
        for i in range(5):
            for t, v in enumerate(C[i]):
                out[(v,)] = out.get((v,), 0) + (1 << (i * h + t))
        return out

    def _prod_poly(self, A, fA, B, fB, tag, depth, toom=0):
        """Return a poly (dict monomial -> int coef) whose value is exactly the
        integer  fA(wv) * fB(wv).  May create auxiliary words + assertions.

        toom > 0                       -> Toom-3 split (5 sub-multiplies of n/3).
        depth <= 0 (or operands short) -> schoolbook bit-pair poly (quadratic).
        depth > 0                      -> a LINEAR poly over freshly materialised
                                          sub-product words (Karatsuba)."""
        na, nb = len(A), len(B)
        if toom > 0 and min(na, nb) > max(self.kmin, 6):
            return self._toom3_poly(A, fA, B, fB, tag, depth, toom - 1)
        if depth <= 0 or min(na, nb) <= self.kmin:
            poly = {}
            for i, u in enumerate(A):
                for j, v in enumerate(B):
                    m = (u, v) if u <= v else (v, u)
                    poly[m] = poly.get(m, 0) + (1 << (i + j))
            return poly

        Q = self.qb
        h = max(na, nb) // 2
        A0, A1 = A[:h], A[h:]
        B0, B1 = B[:h], B[h:]
        msk = (1 << h) - 1
        fA0 = lambda wv, f=fA, m=msk: f(wv) & m
        fA1 = lambda wv, f=fA, s=h: f(wv) >> s
        fB0 = lambda wv, f=fB, m=msk: f(wv) & m
        fB1 = lambda wv, f=fB, s=h: f(wv) >> s
        same = (A is B) or (A == B)

        z0, f0 = self._prod_word(A0, fA0, B0, fB0, tag + 'L', depth - 1, 0)
        z2, f2 = self._prod_word(A1, fA1, B1, fB1, tag + 'H', depth - 1, 0)

        def sumword(lo, hi, flo, fhi, nm):
            n = max(len(lo), len(hi)) + 1
            W = Q.word(nm, n, lambda wv, a=flo, b=fhi: a(wv) + b(wv))
            poly = {}
            for t, v in enumerate(lo):
                poly[(v,)] = poly.get((v,), 0) + (1 << t)
            for t, v in enumerate(hi):
                poly[(v,)] = poly.get((v,), 0) + (1 << t)
            for t, v in enumerate(W):
                poly[(v,)] = poly.get((v,), 0) - (1 << t)
            Q.assert_zero(poly, 0, nm)
            return W, (lambda wv, nm=nm: wv[nm])

        S, fS = sumword(A0, A1, fA0, fA1, f"ks:{tag}")
        if same:
            T, fT = S, fS
        else:
            T, fT = sumword(B0, B1, fB0, fB1, f"kt:{tag}")
        zm, fm = self._prod_word(S, fS, T, fT, tag + 'M', depth - 1, 0)

        # MID = zm - z0 - z2 = A0*B1 + A1*B0  (>= 0, < 2^{na+nb-h})
        nmid = na + nb - h
        nmname = f"kmid:{tag}"
        MID = Q.word(nmname, nmid,
                     lambda wv, a=fm, b=f0, c=f2: a(wv) - b(wv) - c(wv))
        poly = {}
        for t, v in enumerate(zm):
            poly[(v,)] = poly.get((v,), 0) + (1 << t)
        for t, v in enumerate(z0):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        for t, v in enumerate(z2):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        for t, v in enumerate(MID):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        Q.assert_zero(poly, 0, nmname)

        out = {}
        for t, v in enumerate(z0):
            out[(v,)] = out.get((v,), 0) + (1 << t)
        for t, v in enumerate(MID):
            out[(v,)] = out.get((v,), 0) + (1 << (h + t))
        for t, v in enumerate(z2):
            out[(v,)] = out.get((v,), 0) + (1 << (2 * h + t))
        return out

    def _prod_word(self, A, fA, B, fB, tag, depth, toom=0):
        """Materialise the exact integer product as a word (no mod-p reduction)."""
        poly = self._prod_poly(A, fA, B, fB, tag, depth, toom)
        Q = self.qb
        nm = f"kw:{tag}"
        W = Q.word(nm, len(A) + len(B), lambda wv, a=fA, b=fB: a(wv) * b(wv))
        pol = dict(poly)
        for t, v in enumerate(W):
            pol[(v,)] = pol.get((v,), 0) - (1 << t)
        Q.assert_zero(pol, 0, f"kz:{tag}")
        return W, (lambda wv, nm=nm: wv[nm])

    # -------------------------------------------------- multiplication API ---
    def mul_eq(self, tag, A, B, nameA, nameB, terms, const):
        """A*B == const + sum coef*term   (mod p)."""
        if self.kdepth <= 0 and self.toom <= 0:
            return super().mul_eq(tag, A, B, nameA, nameB, terms, const)
        fA = lambda wv, n=nameA: wv[n]
        fB = lambda wv, n=nameB: wv[n]
        poly = dict(self._prod_poly(A, fA, B, fB, self._tag(tag),
                                    self.kdepth, self.toom))
        for coef, bits, _nm in terms:
            for t, v in enumerate(bits):
                poly[(v,)] = poly.get((v,), 0) - coef * (1 << t)

        def val(wv, a=nameA, b=nameB, terms=terms, const=const):
            return wv[a] * wv[b] - const - sum(c * wv[nm] for c, _bits, nm in terms)
        self.congruent(poly, -const, f"eq:{tag}", val)

    # ------------------------------------------------------ generic helpers --
    def lin_word(self, name, terms, const, fn):
        """fresh s-bit word == const + sum coef*term (mod p)."""
        Q = self.qb
        wd = Q.word(name, self.s, fn)
        poly = {}
        for coef, bits, _ in terms:
            for t, v in enumerate(bits):
                poly[(v,)] = poly.get((v,), 0) + coef * (1 << t)
        for t, v in enumerate(wd):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        self.congruent(poly, const, f"lin:{name}",
                       lambda wv, terms=terms, const=const, name=name:
                       const + sum(c * wv[nm] for c, _b, nm in terms) - wv[name])
        return wd

    # ------------------------------------------------ one-hot MUX look-up ----
    def mux_word(self, name, sel, consts, fn):
        """word == consts[t] for the unique t with sel[t] == 1.

        ZERO ancillas.  For every output bit c let  S_c = {t : bit c of consts[t]}.
        Add   sum_{t in S_c} sel_t (1 - w_c)  +  sum_{t not in S_c} sel_t w_c.
        Every summand is a product of two binaries, so the penalty is >= 0, and
        it vanishes exactly when w_c agrees with the selected constant given the
        one-hot constraint sum_t sel_t = 1 (asserted separately).
        Cost: 0 new variables, |sel| couplers per output bit."""
        Q = self.qb
        wd = Q.word(name, self.s, fn)
        for c, v in enumerate(wd):
            for t, u in enumerate(sel):
                m = (u, v) if u < v else (v, u)
                if (consts[t] >> c) & 1:
                    Q._add(Q.pen, (u,), 1)
                    Q._add(Q.pen, m, -1)
                else:
                    Q._add(Q.pen, m, 1)
        return wd

    def sum_word(self, name, sel, selnames, consts, fn):
        """the ORIGINAL look-up: word == sum_t sel_t * consts[t]  (mod p)."""
        Q = self.qb
        wd = Q.word(name, self.s, fn)
        poly = {}
        for t in range(len(sel)):
            poly[(sel[t],)] = poly.get((sel[t],), 0) + consts[t]
        for t, v in enumerate(wd):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        self.congruent(poly, 0, f"lin:{name}",
                       lambda wv, consts=consts, selnames=selnames, name=name:
                       sum(c * wv[nm] for c, nm in zip(consts, selnames)) - wv[name])
        return wd

    def sel_word(self, name, sel, selnames, consts, fn, mux):
        if mux:
            return self.mux_word(name, sel, consts, fn)
        return self.sum_word(name, sel, selnames, consts, fn)

    # ------------------------------------------------- one-hot digit vectors -
    def onehot_square(self, pre, w, digitname):
        """the ORIGINAL selector: D free bits + a squared cardinality constraint.
        D variables, D(D-1)/2 couplers -- the D^2 is what caps the window width."""
        Q = self.qb
        D = 1 << w
        u = [Q.new(f"{pre}_{t}", 'input') for t in range(D)]
        for t, v in enumerate(u):
            Q.trace.append(('word', f"{pre}_{t}", [v],
                            (lambda wv, dn=digitname, t=t: 1 if wv[dn] == t else 0)))
        Q.add_square({v: 1 for v in u}, -1)
        return u

    def _lit_tree(self, pre, w, digitname, lo):
        """one-hot over w bits, built as a tree of ANDs of literals.
        Structurally one-hot -- no cardinality penalty at all."""
        Q = self.qb
        outs = [None]                       # None == the constant 1
        for i in range(w):
            d = Q.new(f"{pre}d{i}", 'input')
            nd = Q.new(f"{pre}n{i}", 'input')
            Q.trace.append(('word', f"{pre}d{i}", [d],
                            (lambda wv, dn=digitname, k=lo + i: (wv[dn] >> k) & 1)))
            Q.trace.append(('word', f"{pre}n{i}", [nd],
                            (lambda wv, dn=digitname, k=lo + i: 1 - ((wv[dn] >> k) & 1)))),
            Q.add_square({d: 1, nd: 1}, -1)          # nd == 1 - d
            zero = [(nd if v is None else Q.AND(v, nd)) for v in outs]
            one = [(d if v is None else Q.AND(v, d)) for v in outs]
            outs = zero + one
        return outs

    def onehot_tree(self, pre, w, digitname):
        """D one-hot indicators from w binary digit bits.

        Split the digit in half and cross-multiply the two literal trees:
        ~D + 4*2^{w/2} AND variables, ~4D couplers, and NO cardinality penalty
        (exactly one indicator is 1 by construction).  This is what makes a wide
        window affordable once the look-up itself is free."""
        Q = self.qb
        if w == 0:
            return [None]
        w1 = w // 2
        w2 = w - w1
        if w1 == 0:
            return self._lit_tree(pre, w, digitname, 0)
        lo = self._lit_tree(pre + "L", w1, digitname, 0)
        hi = self._lit_tree(pre + "H", w2, digitname, w1)
        out = []
        for b in range(1 << w2):
            for a in range(1 << w1):
                out.append(Q.AND(lo[a], hi[b]))
        return out


# ============================================================ windowed comb ==
def build_comb(p, B, table, T, w, chunk=16, mode='binary',
               mux=True, kdepth=0, kmin=8, signed=False, onehot='square',
               verbose=False):
    """Windowed comb, three-multiplication affine step, d != 0 gadget.

    unsigned (signed=False):
        table[j][t] = (t+1) * 2^{wj} G,  t = 0..2^w-1;  T must already be offset
        by  sum_j 2^{wj} G.
    signed (signed=True):
        table[j][t] = (2t+1) * 2^{wj} G, t = 0..2^{w-1}-1, plus one sign bit per
        window.  T must be  2*T_true - (2^{Mw} - 1) G.

    Witness inputs: wv['_u{j}'] = digit index, wv['_s{j}'] = sign bit (signed).
    Returns (L, SEL) with SEL[j] = (one_hot_vars, sign_var_or_None).
    """
    M = len(table)
    D = len(table[0])
    L = Ladder2(p, chunk=chunk, mode=mode, kdepth=kdepth, kmin=kmin)
    Q = L.qb
    s = L.s

    dw = (w - 1) if signed else w
    assert D == 1 << dw
    SEL = []
    for j in range(M):
        if onehot == 'tree':
            u = L.onehot_tree(f"u{j}", dw, f"_u{j}")
        else:
            u = L.onehot_square(f"u{j}", dw, f"_u{j}")
        sg = None
        if signed:
            sg = Q.new(f"sg{j}", 'input')
            Q.trace.append(('word', f"sg{j}", [sg],
                            (lambda wv, j=j: wv[f"_s{j}"])))
        SEL.append((u, sg))

    def lookup(name, j, coord):
        u, _ = SEL[j]
        consts = [table[j][t][coord] % p for t in range(D)]
        names = [f"u{j}_{t}" for t in range(D)]
        fn = (lambda wv, j=j, coord=coord: table[j][wv[f"_u{j}"]][coord] % p)
        return L.sel_word(name, u, names, consts, fn, mux)

    def yword(name, j):
        """window j's y coordinate, sign applied when signed."""
        _, sg = SEL[j]
        ya = lookup(name + "a", j, 1)
        if not signed:
            return ya, name + "a"
        z = [Q.AND(sg, v) for v in ya]
        wd = Q.word(name, s, lambda wv, na=name + "a", j=j:
                    (wv[na] * (1 - 2 * wv[f"_s{j}"])) % p)
        poly = {}
        for t, v in enumerate(ya):
            poly[(v,)] = poly.get((v,), 0) + (1 << t)
        for t, v in enumerate(z):
            poly[(v,)] = poly.get((v,), 0) - 2 * (1 << t)
        for t, v in enumerate(wd):
            poly[(v,)] = poly.get((v,), 0) - (1 << t)
        L.congruent(poly, 0, f"lin:{name}",
                    lambda wv, na=name + "a", j=j, nm=name:
                    wv[na] * (1 - 2 * wv[f"_s{j}"]) - wv[nm])
        return wd, name

    x1 = lookup("x1_0", 0, 0)
    nx = "x1_0"
    y1, ny = yword("y1_0", 0)

    for j in range(1, M):
        x2 = lookup(f"x2_{j}", j, 0)
        n2x = f"x2_{j}"
        y2, n2y = yword(f"y2_{j}", j)
        d = L.lin_word(f"d{j}", [(1, x2, n2x), (-1, x1, nx)], 0,
                       lambda wv, a=n2x, b=nx: (wv[a] - wv[b]) % p)
        e = L.lin_word(f"e{j}", [(1, y2, n2y), (-1, y1, ny)], 0,
                       lambda wv, a=n2y, b=ny: (wv[a] - wv[b]) % p)
        lam = Q.word(f"lam{j}", s, lambda wv, j=j:
                     wv[f"e{j}"] * pow(wv[f"d{j}"], -1, p) % p)
        L.mul_eq(f"lam{j}", lam, d, f"lam{j}", f"d{j}", [(1, e, f"e{j}")], 0)
        for c in (0, p):
            not_equal(L, Q, d, f"d{j}", c, f"ne{j}_{c}")
        x3 = Q.word(f"x3_{j}", s, lambda wv, j=j, nx=nx, n2x=n2x:
                    (wv[f"lam{j}"] ** 2 - wv[nx] - wv[n2x]) % p)
        L.mul_eq(f"x3e{j}", lam, lam, f"lam{j}", f"lam{j}",
                 [(1, x3, f"x3_{j}"), (1, x1, nx), (1, x2, n2x)], 0)
        mm = L.lin_word(f"m{j}", [(1, x1, nx), (-1, x3, f"x3_{j}")], 0,
                        lambda wv, nx=nx, j=j: (wv[nx] - wv[f"x3_{j}"]) % p)
        y3 = Q.word(f"y3_{j}", s, lambda wv, j=j, ny=ny:
                    (wv[f"lam{j}"] * wv[f"m{j}"] - wv[ny]) % p)
        L.mul_eq(f"y3e{j}", lam, mm, f"lam{j}", f"m{j}",
                 [(1, y3, f"y3_{j}"), (1, y1, ny)], 0)
        x1, y1, nx, ny = x3, y3, f"x3_{j}", f"y3_{j}"
        if verbose:
            print(f"  window {j}/{M}: vars={Q.n}", flush=True)

    Tx, Ty = T
    L.congruent({(v,): (1 << t) for t, v in enumerate(x1)}, -Tx, "finalX",
                lambda wv, nx=nx: wv[nx] - Tx)
    L.congruent({(v,): (1 << t) for t, v in enumerate(y1)}, -Ty, "finalY",
                lambda wv, ny=ny: wv[ny] - Ty)
    Q.finalize()
    return L, SEL


# ============================================= x-only (Semaev S_3) chain =====
def s3(X1, X2, X3, Bc, p):
    """third Semaev summation polynomial for y^2 = x^3 + B  (a = 0)."""
    return ((X1 - X2) ** 2 * X3 * X3
            - 2 * ((X1 + X2) * X1 * X2 + 2 * Bc) * X3
            + (X1 * X2) ** 2 - 4 * Bc * (X1 + X2)) % p


def build_semaev(p, Bc, table, xT, w, chunk=16, mode='binary',
                 mux=True, kdepth=0, kmin=8, onehot='square', verbose=False):
    """x-only comb: no y coordinate anywhere, no inversion, no division.

    table[j][t] = x( (2t+1) * 2^{wj} G ),  t = 0..2^{w-1}-1  (magnitudes only;
    x(-P) = x(P), so each digit's sign is free).  Chain of S_3's:

        S_3(z_{j-1}, x_j, z_j) = 0        j = 1..M-2
        S_3(z_{M-2}, x_{M-1}, xT) = 0     final, xT a compile-time constant

    RELAXATION: S_3(a,b,c) = 0 iff c = x(P +- Q), so every step branches two ways.

    Witness inputs: wv['_u{j}'] digit index, wv['_z{j}'] the intermediate
    x-coordinate the chain is supposed to take (j = 1..M-2).
    """
    M = len(table)
    D = len(table[0])
    L = Ladder2(p, chunk=chunk, mode=mode, kdepth=kdepth, kmin=kmin)
    Q = L.qb
    s = L.s

    SEL = []
    for j in range(M):
        if onehot == 'tree':
            u = L.onehot_tree(f"u{j}", w - 1, f"_u{j}")
        else:
            u = L.onehot_square(f"u{j}", w - 1, f"_u{j}")
        SEL.append(u)

    def lookup(name, j):
        consts = [table[j][t] % p for t in range(D)]
        names = [f"u{j}_{t}" for t in range(D)]
        fn = (lambda wv, j=j: table[j][wv[f"_u{j}"]] % p)
        return L.sel_word(name, SEL[j], names, consts, fn, mux)

    zname = "z0"
    z = lookup("z0", 0)

    for j in range(1, M):
        xj = lookup(f"xt{j}", j)
        nxj = f"xt{j}"
        last = (j == M - 1)
        if not last:
            nz = f"z{j}"
            zn = Q.word(nz, s, lambda wv, j=j: wv[f"_z{j}"])
        else:
            nz, zn = None, None

        v = L.lin_word(f"v{j}", [(1, z, zname), (1, xj, nxj)], 0,
                       lambda wv, a=zname, b=nxj: (wv[a] + wv[b]) % p)
        dd = L.lin_word(f"dd{j}", [(1, z, zname), (-1, xj, nxj)], 0,
                        lambda wv, a=zname, b=nxj: (wv[a] - wv[b]) % p)
        u_ = Q.word(f"uu{j}", s, lambda wv, a=zname, b=nxj: wv[a] * wv[b] % p)
        L.mul_eq(f"u{j}", z, xj, zname, nxj, [(1, u_, f"uu{j}")], 0)
        cc = Q.word(f"cc{j}", s,
                    lambda wv, j=j: (wv[f"v{j}"] * wv[f"uu{j}"] + 2 * Bc) % p)
        L.mul_eq(f"c{j}", v, u_, f"v{j}", f"uu{j}", [(1, cc, f"cc{j}")], -2 * Bc)
        uu2 = Q.word(f"u2{j}", s, lambda wv, j=j: wv[f"uu{j}"] ** 2 % p)
        L.mul_eq(f"u2e{j}", u_, u_, f"uu{j}", f"uu{j}", [(1, uu2, f"u2{j}")], 0)

        poly = {}
        if not last:
            g = Q.word(f"g{j}", s, lambda wv, j=j: wv[f"dd{j}"] * wv[f"z{j}"] % p)
            L.mul_eq(f"ge{j}", dd, zn, f"dd{j}", f"z{j}", [(1, g, f"g{j}")], 0)
            gg = Q.word(f"gg{j}", s, lambda wv, j=j: wv[f"g{j}"] ** 2 % p)
            L.mul_eq(f"gge{j}", g, g, f"g{j}", f"g{j}", [(1, gg, f"gg{j}")], 0)
            hh = Q.word(f"hh{j}", s, lambda wv, j=j: wv[f"cc{j}"] * wv[f"z{j}"] % p)
            L.mul_eq(f"hhe{j}", cc, zn, f"cc{j}", f"z{j}", [(1, hh, f"hh{j}")], 0)
            cg, ch = 1, -2
        else:
            gg = Q.word(f"gg{j}", s, lambda wv, j=j: wv[f"dd{j}"] ** 2 % p)
            L.mul_eq(f"gge{j}", dd, dd, f"dd{j}", f"dd{j}", [(1, gg, f"gg{j}")], 0)
            hh = cc
            cg, ch = xT * xT % p, -2 * xT % p
        for t, vv in enumerate(gg):
            poly[(vv,)] = poly.get((vv,), 0) + cg * (1 << t)
        for t, vv in enumerate(hh):
            poly[(vv,)] = poly.get((vv,), 0) + ch * (1 << t)
        for t, vv in enumerate(uu2):
            poly[(vv,)] = poly.get((vv,), 0) + (1 << t)
        for t, vv in enumerate(v):
            poly[(vv,)] = poly.get((vv,), 0) - 4 * Bc * (1 << t)
        nh = (f"hh{j}" if not last else f"cc{j}")
        L.congruent(poly, 0, f"s3:{j}",
                    lambda wv, j=j, Bc=Bc, cg=cg, ch=ch, nh=nh:
                    (cg * wv[f"gg{j}"] + ch * wv[nh]
                     + wv[f"u2{j}"] - 4 * Bc * wv[f"v{j}"]))
        if not last:
            z, zname = zn, nz
        if verbose:
            print(f"  s3 step {j}/{M-1}: vars={Q.n}", flush=True)
    Q.finalize()
    return L, SEL
