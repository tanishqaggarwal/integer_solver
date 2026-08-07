#!/usr/bin/env python3
"""Agent AC. Exact popcount distribution of k uniform on [0,N), secp256k1 N.

Everything here is exact integer / Fraction arithmetic.  No normal approximation,
no Chernoff bound is used to PRODUCE a number -- Chernoff appears only as an
independent bracket in ac_check.py.
"""
from math import comb
from fractions import Fraction
import json, sys

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
P = 2**256 - 2**32 - 977
TWO256 = 2**256

# ---------------------------------------------------------------- digit DP
def count_by_popcount_below(M, nbits=256):
    """Exact vector cnt[b] = #{ k in [0,M) : popcount(k) == b }, via digit DP
    over the binary expansion of M.  O(nbits^2) integer ops."""
    cnt = [0] * (nbits + 1)
    ones_above = 0                      # popcount of M's bits strictly above i
    for i in range(nbits - 1, -1, -1):
        if (M >> i) & 1:
            # k agrees with M above bit i, has 0 where M has 1 at bit i,
            # and is free on the i low bits.  All such k < M.
            for j in range(0, i + 1):
                cnt[ones_above + j] += comb(i, j)
            ones_above += 1
    return cnt

def main():
    cnt = count_by_popcount_below(N)
    tot = sum(cnt)
    assert tot == N, (tot, N)                       # THE verification the brief asked for
    binom = [comb(256, b) for b in range(257)]
    assert sum(binom) == TWO256

    # difference between [0,N) and [0,2^256)
    gap = [binom[b] - cnt[b] for b in range(257)]   # counts of k in [N, 2^256)
    assert sum(gap) == TWO256 - N
    assert all(g >= 0 for g in gap)

    out = {
        "N": str(N),
        "cnt": [str(c) for c in cnt],
        "binom": [str(c) for c in binom],
        "gap": [str(g) for g in gap],
        "sum_ok": tot == N,
        "twofivesix_minus_N": str(TWO256 - N),
    }
    json.dump(out, open("ac_prior.json", "w"))

    print("SUM CHECK: sum_b cnt[b] == N :", tot == N)
    print("2^256 - N =", TWO256 - N, " = 2^%.4f" % (len(bin(TWO256-N))-2 - 1 + 0.0,))
    import math
    print("           log2 =", math.log2(TWO256 - N))
    print()
    print(" b   cnt[b] (exact, [0,N))            C(256,b)                       rel.deficit")
    for b in list(range(0, 20)) + [24, 56, 100, 126, 127, 128, 129, 130, 192, 240, 246, 250, 255, 256]:
        d = Fraction(gap[b], binom[b]) if binom[b] else 0
        print("%4d %-30d %-30d %.6e" % (b, cnt[b], binom[b], float(d)))
    print()
    mx = max(range(257), key=lambda b: Fraction(gap[b], binom[b]))
    print("largest relative deficit over all b: b=%d  gap/binom = %.6e = 2^%.3f"
          % (mx, float(Fraction(gap[mx], binom[mx])),
             math.log2(float(Fraction(gap[mx], binom[mx])))))
    # weighted: total deficit fraction
    import math
    print("total deficit (2^256-N)/2^256 = 2^%.4f" % math.log2((TWO256-N)/TWO256))

if __name__ == "__main__":
    main()
