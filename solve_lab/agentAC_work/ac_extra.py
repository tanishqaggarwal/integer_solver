#!/usr/bin/env python3
"""Agent AC. Supplementary exact numbers: the full P(w<=B) curve at the B values
that matter operationally, the bracket width once families C-H are included,
and an mpmath cross-check of every headline exponent."""
from math import comb
import math, json
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 60

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
TWO256 = 2**256

def count_by_popcount_below(M, nbits=256):
    cnt = [0]*(nbits+1); oa = 0
    for i in range(nbits-1, -1, -1):
        if (M >> i) & 1:
            for j in range(i+1): cnt[oa+j] += comb(i, j)
            oa += 1
    return cnt

cnt = count_by_popcount_below(N); assert sum(cnt) == N
binom = [comb(256,b) for b in range(257)]
rem = [int(x) for x in json.load(open("ac_cond.json"))["rem"]]
LO, HI = 10, 246
Z = sum(cnt[LO:HI+1])

def lg_dec(num, den):
    """log2 of num/den to ~40 significant digits, independent of math.log."""
    # use Decimal ln
    return float((Decimal(num).ln() - Decimal(den).ln()) / Decimal(2).ln())

def lg(fr):
    return math.log(fr.numerator, 2) - math.log(fr.denominator, 2)

def cum1(B):
    if B < LO: return Fraction(0)
    return Fraction(sum(cnt[LO:min(B,HI)+1]), Z)

print("EXACT POSTERIOR TAIL  P(w <= B | 10<=w<=246), and the rigorous bracket once")
print("families C-H (signed<=7, BSGS, window, small-multiple, lambda-box, orbit) are")
print("also removed.  Bracket width is reported in bits of the exponent.")
print()
print("  B   |  meaning                                  |  P(w<=B)        | bracket width (bits)")
labels = {14:"reach of THIS BOX (AB S3: 2^47 time, 2^30 mem)",
          20:"", 24:"'actionable' band (AB S5)",
          32:"", 40:"", 48:"",
          52:"memory-aware rho crossover at 2^30 memory",
          56:"rho crossover at 2^40 memory (AB S5 payoff band)",
          64:"rho crossover at 2^60 memory",
          80:"", 104:"AB Thm B unbounded-memory crossover (struck)",
          128:"the mode", 148:"AB Thm B break-even B"}
for B in sorted(labels):
    p = cum1(B)
    num_lo = sum(max(0,cnt[b]-rem[b]) for b in range(LO,B+1))
    den_lo = sum(max(0,cnt[b]-rem[b]) for b in range(LO,HI+1))
    plo = Fraction(num_lo, Z); phi = Fraction(sum(cnt[LO:B+1]), den_lo)
    w = lg(phi) - lg(plo)
    print("%5d | %-41s | 2^%-11.4f | %.2e" % (B, labels[B], lg(p), w))

print()
print("independent Decimal recomputation of the three headline exponents:")
for B in [14,24,56]:
    num = sum(cnt[LO:B+1]); den = Z
    print("   B=%3d  math.log route 2^%.6f   Decimal route 2^%.6f   agree to %.1e"
          % (B, lg(Fraction(num,den)), lg_dec(num,den), abs(lg(Fraction(num,den))-lg_dec(num,den))))

print()
print("complementary framing (the number the user actually asked for):")
for B in [14,24,56]:
    p = cum1(B)
    print("   P(w > %3d) = 1 - 2^%.4f = 1 - %.3e   i.e. w is ABOVE the band with"
          " certainty 1 - 2^%.2f" % (B, lg(p), float(p), lg(p)))

print()
print("expected number of instances you would have to draw before seeing w <= B:")
for B in [14,24,56]:
    p = cum1(B)
    print("   B=%3d : 2^%.2f draws  (age of the universe in Planck times ~ 2^202)" % (B, -lg(p)))
