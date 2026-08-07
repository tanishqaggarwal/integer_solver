#!/usr/bin/env python3
"""Agent AC. Exact prior, exact posterior, quantile table, headline tails,
and a second-route bracket on every tail so no exponent rests on one sum."""
from math import comb
import math, json
from fractions import Fraction

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
TWO256 = 2**256

def count_by_popcount_below(M, nbits=256):
    cnt = [0] * (nbits + 1)
    ones_above = 0
    for i in range(nbits - 1, -1, -1):
        if (M >> i) & 1:
            for j in range(0, i + 1):
                cnt[ones_above + j] += comb(i, j)
            ones_above += 1
    return cnt

def lg(fr):
    """log2 of a positive Fraction, exactly enough (CPython math.log is exact on big ints)."""
    if fr == 0:
        return float("-inf")
    return math.log(fr.numerator, 2) - math.log(fr.denominator, 2)

cnt = count_by_popcount_below(N)
assert sum(cnt) == N
binom = [comb(256, b) for b in range(257)]

LO, HI = 10, 246                      # w >= 10 (X), w <= 246 (Y)
Z = sum(cnt[LO:HI + 1])               # exact posterior normaliser

def P0(b):  return Fraction(cnt[b], N)
def P1(b):  return Fraction(cnt[b], Z) if LO <= b <= HI else Fraction(0)

def cum0(B): return Fraction(sum(cnt[:B + 1]), N)
def cum1(B): return Fraction(sum(cnt[max(LO,0):min(B, HI) + 1]), Z) if B >= LO else Fraction(0)

print("=" * 78)
print("1. EXACT PRIOR  (k uniform on [0,N))")
print("=" * 78)
print("sum_b cnt[b] == N :", sum(cnt) == N)
gap = [binom[b] - cnt[b] for b in range(257)]
first = next(b for b in range(257) if gap[b] != 0)
print("first popcount class where [0,N) differs from [0,2^256):  b =", first,
      "  (gap =", gap[first], ")")
print("=> for every b <= %d, #{k<N : popcount=b} == C(256,b) EXACTLY (deficit 0, not small)" % (first - 1))
print("total deficit (2^256-N)/2^256 = 2^%.4f" % lg(Fraction(TWO256 - N, TWO256)))
worst = max(range(LO, HI + 1), key=lambda b: Fraction(gap[b], binom[b]))
print("worst relative deficit inside the live range [10,246]: b=%d, %.4e = 2^%.2f"
      % (worst, float(Fraction(gap[worst], binom[worst])), lg(Fraction(gap[worst], binom[worst]))))

print()
print("=" * 78)
print("2. POSTERIOR after w>=10 (X) and w<=246 (Y)")
print("=" * 78)
excl = Fraction(N - Z, N)
print("P_prior(excluded) = P(w<=9) + P(w>=247) = %s" % excl)
print("                  = 2^%.4f" % lg(excl))
print("TV(prior, posterior) = P_prior(excluded) = 2^%.4f" % lg(excl))
print("information gained   = -log2(1-P_excl)  = %.4e bits = 2^%.2f bits"
      % (-math.log2(1 - float(excl)) if float(excl) else 0.0, lg(excl / Fraction(math.log(2)).limit_denominator(10**12))))
kl = -float(excl) / math.log(2)   # ~ log2(1/(1-eps))
print("P(w<=9)   = 2^%.4f      P(w>=247) = 2^%.4f" % (lg(cum0(9)), lg(Fraction(sum(cnt[247:]), N))))
print("(these are NOT equal: [0,N) is not symmetric.  Under [0,2^256) they would be.)")
sym = Fraction(sum(binom[:10]), TWO256), Fraction(sum(binom[247:]), TWO256)
print("   under the untruncated model both are 2^%.4f" % lg(sym[0]))

# --- cross-check against agent X's per-bit 2^-200.8 ---------------------------
fam = {
  "weight<=9 (X)":            sum(binom[:10]),
  "complement weight<=9 (Y)": sum(binom[247:]),
  "BSGS both ends 2*2^52":    2 * 2**52,
  "a+b*lambda box":           (2 * 2**21) ** 2,
  "34-bit window":            223 * 2**33,
  "small multiple m<=1e7":    256 * 10**7,
  "endomorphism orbit":       1536,
  "NAF weight <= 7 (set)":    733018571531264,
}
U = sum(fam.values())
print()
print("union bound over ALL excluded families:")
for k, v in sorted(fam.items(), key=lambda kv: -kv[1]):
    print("   %-28s %-22d 2^%.2f" % (k, v, math.log2(v)))
print("   %-28s %-22d 2^%.2f" % ("UNION (upper bound)", U, math.log2(U)))
print("fraction of a single-bit fibre (2^255): 2^%.4f    [X reports 2^-200.80]"
      % (math.log2(U) - 255))
print("fraction of the whole space  (2^256):  2^%.4f" % (math.log2(U) - 256))
Ux = U - fam["complement weight<=9 (Y)"] - fam["NAF weight <= 7 (set)"]
print("X's own family list (no Y complement, X's overcounted signed size): 2^%.4f of 2^255"
      % (math.log2(Ux + 1708592205206529) - 255))

print()
print("=" * 78)
print("3. QUANTILES")
print("=" * 78)
eps_list = [("1e-1", Fraction(1,10)), ("1e-2", Fraction(1,100)), ("1e-3", Fraction(1,1000)),
            ("1e-6", Fraction(1,10**6)), ("2^-20", Fraction(1,2**20)),
            ("2^-40", Fraction(1,2**40)), ("2^-80", Fraction(1,2**80))]

def smallest_B(cum, target):
    for B in range(0, 257):
        if cum(B) >= target:
            return B
    return None

def two_sided(cum, eps):
    half = eps / 2
    # L = smallest L with P(w < L) <= half  -> largest L such that cum(L-1) <= half
    L = 0
    for b in range(0, 257):
        if cum(b - 1 if b else -1) if False else cum(b) > half:
            L = b
            break
    # recompute properly: L is the smallest b with cum(b-1) <= half < cum(b) ... we want
    # P(w<L) <= half, maximal L
    L = 0
    for b in range(0, 257):
        c = cum(b - 1) if b > 0 else Fraction(0)
        if c <= half:
            L = b
        else:
            break
    U_ = 256
    for b in range(256, -1, -1):
        tail = Fraction(1) - cum(b)
        if tail <= half:
            U_ = b
        else:
            break
    return L, U_

def cum0m1(B):
    return cum0(B) if B >= 0 else Fraction(0)
def cum1m1(B):
    return cum1(B) if B >= 0 else Fraction(0)

print("PRIOR (uniform on [0,N)) :")
print(" eps      one-sided B: P(w<=B)>=1-eps    two-sided [L,U] with each tail <= eps/2")
for nm, e in eps_list:
    B = smallest_B(cum0m1, 1 - e)
    L, U_ = two_sided(cum0m1, e)
    print("  %-6s  B = %3d   (P(w<=B)-(1-eps) ok)      [%3d, %3d]  width %d" % (nm, B, L, U_, U_ - L))
print()
print("POSTERIOR (conditioned on 10 <= w <= 246) :")
for nm, e in eps_list:
    B = smallest_B(cum1m1, 1 - e)
    L, U_ = two_sided(cum1m1, e)
    print("  %-6s  B = %3d                              [%3d, %3d]  width %d" % (nm, B, L, U_, U_ - L))

print()
print("=" * 78)
print("4. THE HEADLINE TAILS  (AB's cost bands)")
print("=" * 78)
def bracket(B):
    """Independent bracket on S = sum_{b<=B} C(256,b), B < 128:
       C(256,B) <= S <= C(256,B)/(1-r), r = B/(256-B+1)."""
    r = Fraction(B, 256 - B + 1)
    S = sum(binom[:B + 1])
    lo, hi = binom[B], Fraction(binom[B], 1 - r)
    return S, lo, hi, (lo <= S <= hi)

print(" B     P_prior(w<=B)      P_post(w<=B)     ratio-bracket   Chernoff(2^-256*2^(256 H(B/256))) ")
for B in [14, 24, 56, 64, 104]:
    p0, p1 = cum0(B), cum1(B)
    S, lo, hi, ok = bracket(B)
    x = B / 256.0
    H = -(x * math.log2(x) + (1 - x) * math.log2(1 - x))
    ch = 256 * (H - 1)          # log2 upper bound on P(Bin(256,1/2)<=B)
    print("%4d   2^%-12.4f   2^%-12.4f   %s   <= 2^%.4f  (exact 2^%.4f)"
          % (B, lg(p0), lg(p1), "OK" if ok else "FAIL", ch, lg(Fraction(S, TWO256))))

print()
for B in [14, 24, 56]:
    p1 = cum1(B)
    print("P(w <= %3d | 10<=w<=246) = %s" % (B, p1))
    print("      = 2^%.4f  =  %.6e" % (lg(p1), float(p1)))

print()
print("odds form:")
for B in [14, 24, 56]:
    p1 = cum1(B)
    print("   1 in 2^%.2f" % (-lg(p1)))

json.dump({"cnt": [str(c) for c in cnt], "Z": str(Z),
           "P_le": {str(B): [str(cum0(B).numerator), str(cum0(B).denominator)] for B in range(257)}},
          open("ac_post.json", "w"))
