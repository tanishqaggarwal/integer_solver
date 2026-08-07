#!/usr/bin/env python3
"""Agent AC. Exact per-popcount-class sizes of every family the campaign has
excluded, so that "how much did the campaign move the posterior on w" is a
computed number and not an assertion.

Families (from agentX_work/K_CONSTRAINTS.md §1, each re-derived here as a SET;
I take the *exhaustion* claims on the record's authority and compute only the
combinatorics):

  A  unsigned weight <= 9                       (X, exhaustive)      -> kills classes b<=9
  B  complement weight <= 9  i.e. w >= 247      (Y, exhaustive)      -> kills classes b>=247
  C  signed-digit (NAF) weight <= 7             (X, exhaustive)
  D  k < 2^52  or  N-k < 2^52                   (X, BSGS)
  E  ON-bits inside some 34-bit window          (Q, audited)
  F  k = 2^i/m mod N, m <= 10^7                 (Q, audited)
  G  k = a + b*lambda, |a|,|b| < 2^21           (Q, audited)
  H  endomorphism orbit, 1536 points            (Q/X)

A and B are the only ones that are *events about w itself*; C..H are events
about k that cut across weight classes.  Nothing here is multiplied: the
posterior is the prior restricted to the complement of the union, so a UNION
BOUND on the removed count per class is a rigorous bound on the movement.
"""
from math import comb
import json, sys

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
TWO256 = 2**256

# ------------------------------------------------------------------ helpers
def count_by_popcount_below(M, nbits=256):
    cnt = [0] * (nbits + 1)
    ones_above = 0
    for i in range(nbits - 1, -1, -1):
        if (M >> i) & 1:
            for j in range(0, i + 1):
                cnt[ones_above + j] += comb(i, j)
            ones_above += 1
    return cnt

def nafweight(v):
    """NAF weight via the standard identity: nonzero NAF digits of v sit exactly
    at the positions where 3v and v differ."""
    return bin((3 * v) ^ v).count("1")

def nafweight_ref(v):
    """Independent reference implementation (Reitwiesner's algorithm)."""
    w = 0
    while v:
        if v & 1:
            z = 2 - (v & 3)
            v -= z
            w += 1
        v >>= 1
    return w

# ------------------------------------------------- family C: NAF weight <= 7
def naf_popcount_table(nbits, maxnaf):
    """Exact joint census over k in [0, 2^nbits):
         tab[b][m] = #{k : popcount(k)=b, nafweight(k)=m}   (m capped, see below)
    DP over bits LSB->MSB carrying:
       c3   : carry of the addition v + 2v = 3v
       prev : previous bit of v (needed because 2v's bit i is v's bit i-1)
    A NAF digit is nonzero at position i iff bit_i(3v) != bit_i(v).
    Only m <= maxnaf is tracked; anything larger is dropped (we only need the
    small-m classes)."""
    # state: (c3, prev) -> {(b, m): count}
    from collections import defaultdict
    cur = defaultdict(int)
    cur[(0, 0)] = 1
    states = {(0, 0): {(0, 0): 1}}
    st = {(0, 0): {(0, 0): 1}}
    for i in range(nbits):
        nxt = {}
        for (c3, prev), d in st.items():
            for bit in (0, 1):
                # bit i of 3v = bit i of (v + 2v):  v_i + v_{i-1} + c3
                s = bit + prev + c3
                b3 = s & 1
                c3n = s >> 1
                nz = 1 if b3 != bit else 0
                key = (c3n, bit)
                tgt = nxt.setdefault(key, {})
                for (b, m), cnt in d.items():
                    m2 = m + nz
                    if m2 > maxnaf:
                        continue
                    k2 = (b + bit, m2)
                    tgt[k2] = tgt.get(k2, 0) + cnt
        st = nxt
    # flush the remaining carry: after the last bit, 3v may still have high bits
    # (positions nbits, nbits+1) contributed by prev and c3.  Those positions
    # have v-bit 0, so a nonzero 3v bit there is a nonzero NAF digit.
    out = {}
    for (c3, prev), d in st.items():
        tail = prev + c3          # value of the remaining high part of 3v
        extra = bin(tail).count("1")
        for (b, m), cnt in d.items():
            m2 = m + extra
            if m2 > maxnaf:
                continue
            out[(b, m2)] = out.get((b, m2), 0) + cnt
    tab = [[0] * (maxnaf + 1) for _ in range(nbits + 1)]
    for (b, m), c in out.items():
        tab[b][m] += c
    return tab

def selftest_naf():
    # identity check
    for v in range(0, 5000):
        assert nafweight(v) == nafweight_ref(v), v
    for v in [N, N - 1, TWO256 - 1, 2**255 + 7, 12345678901234567890]:
        assert nafweight(v) == nafweight_ref(v), v
    # DP check by brute force on 14 bits
    nb = 14
    tab = naf_popcount_table(nb, 6)
    brute = [[0] * 7 for _ in range(nb + 1)]
    for v in range(2**nb):
        m = nafweight(v)
        if m <= 6:
            brute[bin(v).count("1")][m] += 1
    assert tab == brute, "NAF DP disagrees with brute force"
    print("selftest: NAF identity ok (5000 values + 5 big); DP == brute force on 2^14  OK")

# ---------------------------------------------------------------------- main
def main():
    selftest_naf()
    cntN = count_by_popcount_below(N)
    assert sum(cntN) == N
    print()

    # ---- NAF weight of N itself (needed for the wrap argument) -------------
    wN = nafweight(N)
    assert wN == nafweight_ref(N)
    print("popcount(N) =", bin(N).count("1"), "   NAFweight(N) =", wN,
          "  (Reitwiesner minimum signed weight of N)")

    # ---- family C ----------------------------------------------------------
    tab = naf_popcount_table(256, 7)
    C_b = [sum(tab[b][:8]) for b in range(257)]     # NAF weight <= 7
    C_tot = sum(C_b)
    enum_size = sum(comb(256, m) * 2**m for m in range(0, 8))
    print("family C  |{k<2^256 : NAFwt(k)<=7}| = %d = 2^%.2f" % (C_tot, __import__('math').log2(C_tot)))
    print("          (X's enumeration size sum_{m<=7} C(256,m)2^m = %d = 2^%.2f -- an overcount,"
          % (enum_size, __import__('math').log2(enum_size)))
    print("           since representations are not unique; the SET is what matters)")

    # ---- family D ----------------------------------------------------------
    D_lo = [comb(52, b) if b <= 52 else 0 for b in range(257)]
    hi = count_by_popcount_below(N)
    lo = count_by_popcount_below(N - 2**52 + 1)
    D_hi = [hi[b] - lo[b] for b in range(257)]      # k in [N-2^52+1, N)
    assert all(x >= 0 for x in D_hi)
    assert sum(D_hi) == 2**52 - 1

    # ---- family E: bits inside a 34-bit window -----------------------------
    # k = a*2^s, a odd, a < 2^34, s in [0,222]  (unique representation)
    E_b = [223 * comb(33, b - 1) if 1 <= b <= 34 else 0 for b in range(257)]

    # ---- families F,G,H: no useful per-class structure; use their total size
    F_tot = 256 * 10**7
    G_tot = (2 * 2**21) ** 2
    H_tot = 1536
    FGH = F_tot + G_tot + H_tot

    # ---- assemble a per-class UPPER BOUND on the removed count -------------
    rem = [0] * 257
    for b in range(257):
        r = C_b[b] + D_lo[b] + D_hi[b] + E_b[b] + FGH
        rem[b] = min(r, cntN[b])
    # classes b<=9 and b>=247 are removed entirely by A and B
    for b in range(0, 10):
        rem[b] = cntN[b]
    for b in range(247, 257):
        rem[b] = cntN[b]

    import math
    print()
    print("  b     cnt[b]            removed<=          rel   log2(rel)")
    for b in [10, 11, 12, 14, 16, 20, 24, 30, 36, 40, 48, 56, 64, 80, 100, 128,
              160, 192, 200, 210, 220, 230, 240, 246]:
        rel = rem[b] / cntN[b]
        print("%4d  %-18.6e %-18.6e %.3e  %8.2f" %
              (b, float(cntN[b]), float(rem[b]), rel, math.log2(rel) if rel else float('-inf')))

    # the wrap argument for family C (see W_POSTERIOR.md S2.3)
    print()
    print("wrap bound: any k = N - u with NAFwt(u)<=7 and popcount(k)=b gives a signed")
    print("            representation of N with <= b+7 nonzero digits, so b >= %d - 7 = %d."
          % (wN, wN - 7))

    json.dump({"cntN": [str(x) for x in cntN],
               "C_b": [str(x) for x in C_b],
               "D_lo": [str(x) for x in D_lo],
               "D_hi": [str(x) for x in D_hi],
               "E_b": [str(x) for x in E_b],
               "FGH": str(FGH),
               "rem": [str(x) for x in rem],
               "nafN": wN, "popcountN": bin(N).count("1")},
              open("ac_cond.json", "w"))

if __name__ == "__main__":
    main()
