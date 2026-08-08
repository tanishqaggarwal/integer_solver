#!/usr/bin/env python3
"""degree.py -- how the quadratic form is reached: AND ancillas vs Ishikawa.

Two measurements:

  (1) AND CACHE HIT RATE.  Every partial product a_i b_j is linearised once and
      reused.  How much reuse is actually there?  For a general product: none
      (all s^2 pairs are distinct).  For a squaring: exactly half, plus the
      diagonal, which is why a*a costs s(s+1)/2 - s ANDs.  Measured, per
      variant, at s = 256.

  (2) ISHIKAWA.  The alternative is not to linearise the products at all and
      instead square the column form directly, which is a degree-4 polynomial,
      then reduce each high-order monomial by Ishikawa's method (1 ancilla per
      quartic monomial, 1 per cubic).  We count those monomials exactly, on the
      same column structure, without building the (astronomically large) form.
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmqb import MMQB                                            # noqa: E402
from mm import build_modmul                                      # noqa: E402
from measure import P, S                                         # noqa: E402


class Dry(MMQB):
    """records the column bags of every identity instead of compiling them."""

    def __init__(self, **kw):
        super().__init__(**kw)
        self.bags = []

    def assert_terms(self, terms, consts, tag):
        from mmqb import best_split
        acc = defaultdict(int)
        for mono, sg, sh in terms:
            acc[tuple(sorted(set(mono)))] += sg << sh
        pos, neg = defaultdict(list), defaultdict(list)
        for mono, coef in acc.items():
            if not coef:
                continue
            for sg, sh in best_split(coef):
                (pos if sg > 0 else neg)[sh].append(mono)
        self.bags.append((tag, pos, neg))

    def mono_var(self, mono):            # never called in dry mode
        raise RuntimeError


def columns(mult='schoolbook', leaf=24, red='naf', square=False, s=S, p=P):
    Q = Dry(mode='wallace')
    A = Q.mkword('A', s, lambda wv: 0)
    B = A if square else Q.mkword('B', s, lambda wv: 0)
    C = Q.mkword('C', s, lambda wv: 0)
    build_modmul(Q, p, A, B, C, mult=mult, leaf=leaf, red=red)
    return Q


def ishikawa_count(Q):
    """monomials of degree 3 and 4 produced by squaring each column form.

    A column of the identity is balanced as  (sum of terms + carries - 2*cout)^2.
    Squaring a form whose atoms are the degree-2 monomials a_i b_j produces one
    cross term per unordered pair of atoms; a pair of distinct products in the
    same column shares no variable (i+j = k+l and i = k forces j = l), so every
    cross term is a distinct quartic.  Cross terms between a product and a
    linear atom (word bit, quotient bit, carry) are cubic."""
    quartic = cubic = 0
    for tag, pos, neg in Q.bags:
        cols = defaultdict(list)
        for d in (pos, neg):
            for c, l in d.items():
                cols[c] += l
        # carry width per column, same interval rule as the binary balancer
        cin = 0
        for c in sorted(cols):
            atoms = cols[c]
            deg2 = [m for m in atoms if len(m) == 2]
            deg1 = [m for m in atoms if len(m) == 1]
            n2, n1 = len(deg2), len(deg1) + cin
            nb = max(0, ((len(atoms) + cin) // 2).bit_length())
            n1 += nb
            quartic += n2 * (n2 - 1) // 2
            cubic += n2 * n1
            cin = nb
    return quartic, cubic


if __name__ == '__main__':
    print("=" * 86)
    print("(1) AND CACHE HIT RATE at s = 256")
    print("=" * 86)
    print(f"{'variant':>40} {'lookups':>10} {'hits':>10} {'rate':>7} {'AND vars':>10}")
    from measure import modmul
    for key, kw in (("general, schoolbook", dict(mult='schoolbook')),
                    ("general, karatsuba(leaf=24)", dict(mult='karatsuba', leaf=24)),
                    ("SQUARING, schoolbook", dict(mult='schoolbook', square=True)),
                    ("SQUARING, karatsuba(leaf=24)", dict(mult='karatsuba', leaf=24,
                                                          square=True))):
        r = modmul(red='naf', mode='wallace', **kw)
        print(f"{key:>40} {r['and_lookups']:10,d} {r['and_hits']:10,d} "
              f"{100*r['and_hits']/max(1,r['and_lookups']):6.1f}% {r['ands']:10,d}")
    print("""
    A general 256x256 product has s^2 distinct index pairs and no reuse at all;
    the cache's only job there is to keep a repeated monomial in the SAME
    equation from being linearised twice.  A squaring hits on exactly the
    transposed pairs: s(s+1)/2 distinct products, s of them diagonal and free,
    so 32,640 ANDs instead of 65,536 -- already exploited, not an opportunity.
    Karatsuba on a squaring is better still: every sub-product in the tree is
    itself a squaring (P0 = a0^2, P1 = a1^2, P2 = (a0+a1)^2), so the AND count
    falls to 10,678.""")
    print()
    print("=" * 86)
    print("(2) ISHIKAWA vs ONE AND ANCILLA PER PRODUCT, s = 256")
    print("=" * 86)
    print(f"{'variant':>34} {'AND ancillas':>14} {'Ishikawa quartic':>18} "
          f"{'cubic':>12} {'Ishikawa total':>16}")
    for key, kw in (("schoolbook", dict(mult='schoolbook')),
                    ("karatsuba(leaf=24)", dict(mult='karatsuba', leaf=24)),
                    ("schoolbook, squaring", dict(mult='schoolbook', square=True))):
        Q = columns(red='naf', **kw)
        q4, q3 = ishikawa_count(Q)
        nand = len({m for _t, pos, neg in Q.bags for d in (pos, neg)
                    for l in d.values() for m in l if len(m) == 2})
        print(f"{key:>34} {nand:14,d} {q4:18,d} {q3:12,d} {q4 + q3:16,d}")
    print("""
    Ishikawa's reduction spends one ancilla per high-order monomial, and the
    monomials it would have to reduce are the cross terms of the column squares
    -- one per unordered pair of partial products in a column.  That is
    quadratic in the column height where the AND ancilla is linear, and it is
    also the wrong shape for the hardware: the reduced form still contains one
    clique per column.  The AND-ancilla-per-product route is not a default here,
    it is ~85x smaller.

    There is also no better factorisation of the partial-product set to find.
    Sharing requires two partial products to be the same monomial; in a general
    product all s^2 are distinct by construction, so the only structural sharing
    available is the transpose symmetry of a squaring (already taken) and the
    reuse of whole sub-products, which is exactly what Karatsuba/Toom do.""")
