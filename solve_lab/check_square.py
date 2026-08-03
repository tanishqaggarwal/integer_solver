#!/usr/bin/env python3
"""Are the big verifier-check atoms (41285, 41627, 42093, 44941, 30693, 40782,...)
perfect squares Q^2 of a simpler polynomial Q? If so, the check reduces to Q=0 (half
the degree) -- a real simplification. Test by: (1) is every monomial a product of two
'base' monomials? (2) reconstruct Q from the square-root of diagonal terms and verify
Q^2 == atom exactly."""
import json
from collections import defaultdict
from math import isqrt, gcd
from propagate import load_atoms, atom_vars

def poly_mul(a, b):
    out = defaultdict(int)
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(sorted(m1 + m2))
            out[m] += c1 * c2
    return {m: c for m, c in out.items() if c}

def try_sqrt(poly):
    """Attempt to write poly = Q^2. Heuristic: base monomials = sqrt of perfect-square
    diagonal monomials (m where m = d+d). Build Q greedily, verify."""
    # find candidate 'squared' monomials: those that are a doubling of a smaller monomial
    dieg = {}
    for m, c in poly.items():
        if len(m) % 2 == 0:
            half = m[::2]
            if tuple(sorted(half + half)) == tuple(sorted(m)) and c > 0:
                s = isqrt(c)
                if s * s == c:
                    dieg[half] = s
    if not dieg:
        return None
    # Q = sum s_i * half_i  (signs unknown; fix first positive, solve others via cross terms)
    bases = list(dieg)
    # set up signs: pick sign of each base by matching a cross term with base[0]
    sign = {bases[0]: 1}
    for b in bases[1:]:
        cross = tuple(sorted(bases[0] + b))
        cc = poly.get(cross, 0)
        expected = 2 * dieg[bases[0]] * dieg[b]
        if expected == 0:
            sign[b] = 1
        else:
            if cc == expected: sign[b] = 1
            elif cc == -expected: sign[b] = -1
            else: return None
    Q = {b: sign[b] * dieg[b] for b in bases}
    if poly_mul(Q, Q) == poly:
        return Q
    return None

def main():
    A = load_atoms()
    for a in (41285, 41627, 42093, 44941, 30693, 40782, 22049):
        poly = A[a]
        deg = max(len(m) for m in poly)
        Q = try_sqrt(poly)
        if Q:
            print(f"atom {a}: PERFECT SQUARE! Q has {len(Q)} terms (deg {max(len(m) for m in Q)}), vars {sorted(atom_vars(Q))[:8]}")
        else:
            # partial info: how many monomials are 'squared' diagonal
            ndiag = sum(1 for m in poly if len(m)%2==0 and tuple(sorted(m[::2]+m[::2]))==tuple(sorted(m)))
            print(f"atom {a}: not a clean square (deg {deg}, {len(poly)} terms, {ndiag} diagonal-square monomials)")
    # Also: is any equation a square of a KNOWN atom combination? check gcd of coeffs
    print("\ncoefficient gcd of big checks:")
    for a in (41285, 41627, 42093, 44941):
        g = 0
        for c in A[a].values(): g = gcd(g, abs(c))
        print(f"  atom {a}: gcd={g}")

if __name__ == '__main__':
    main()
