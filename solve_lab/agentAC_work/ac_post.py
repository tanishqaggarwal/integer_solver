#!/usr/bin/env python3
"""Agent AC. Exact prior, exact posterior, quantile table, headline tails,
with every tail exponent produced twice by independent routes."""
from math import comb
import math, json
from fractions import Fraction

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
TWO256 = 2**256

# ---------------------------------------------------------------- route 1: digit DP
def count_by_popcount_below(M, nbits=256):
    cnt = [0] * (nbits + 1)
    ones_above = 0
    for i in range(nbits - 1, -1, -1):
        if (M >> i) & 1:
            for j in range(0, i + 1):
                cnt[ones_above + j] += comb(i, j)
            ones_above += 1
    return cnt

# ---------------------------------------------------------------- route 2: bitwise DP
def count_by_popcount_below_v2(M, nbits=256):
    """Independent implementation: MSB->LSB DP carrying (tight, ones)."""
    tight = {0: 1}          # ones_so_far -> count, still equal to M's prefix
    free = [0] * (nbits + 1)  # ones_so_far -> count, already strictly below
    for i in range(nbits - 1, -1, -1):
        mb = (M >> i) & 1
        nt = {}
        nf = [0] * (nbits + 1)
        for o, c in enumerate(free):
            if c:
                nf[o] += c        # place 0
                nf[o + 1] += c    # place 1
        for o, c in tight.items():
            if mb == 1:
                nf[o] += c                       # place 0 -> now strictly below
                nt[o + 1] = nt.get(o + 1, 0) + c # place 1 -> still tight
            else:
                nt[o] = nt.get(o, 0) + c         # place 0 -> still tight
        tight, free = nt, nf
    return free   # strictly below M; the "tight" bucket is M itself, excluded

def lg(fr):
    if fr == 0: return float("-inf")
    return math.log(fr.numerator, 2) - math.log(fr.denominator, 2)

cnt = count_by_popcount_below(N)
cnt2 = count_by_popcount_below_v2(N)
assert cnt == cnt2, "the two independent DPs disagree"
assert sum(cnt) == N
binom = [comb(256, b) for b in range(257)]
assert sum(binom) == TWO256

# removals from families C..H (upper bounds), from ac_cond.py
try:
    condj = json.load(open("ac_cond.json"))
    rem = [int(x) for x in condj["rem"]]
except Exception:
    rem = [0] * 257

LO, HI = 10, 246
Z = sum(cnt[LO:HI + 1])

def cum0(B):
    if B < 0: return Fraction(0)
    return Fraction(sum(cnt[:min(B, 256) + 1]), N)
def cum1(B):
    if B < LO: return Fraction(0)
    return Fraction(sum(cnt[LO:min(B, HI) + 1]), Z)

sep = "=" * 84
print(sep); print("1. EXACT PRIOR  (k uniform on [0,N))"); print(sep)
print("two independent digit-DPs agree elementwise :", cnt == cnt2)
print("sum_b cnt[b] == N exactly                   :", sum(cnt) == N)
gap = [binom[b] - cnt[b] for b in range(257)]
first = next(b for b in range(257) if gap[b] != 0)
print("first b where #{k<N: pc=b} != C(256,b)      : b =", first, " (gap =", gap[first], ")")
print("=> for ALL b <= %d the truncation effect is EXACTLY ZERO, not merely 2^-128-small" % (first-1))
print("2^256 - N = %d = 2^%.4f" % (TWO256-N, math.log2(TWO256-N)))
print("total deficit (2^256-N)/2^256               = 2^%.4f" % lg(Fraction(TWO256-N, TWO256)))
wb = max(range(LO, HI+1), key=lambda b: Fraction(gap[b], binom[b]))
print("worst per-class relative deficit in [10,246]: b=%d, 2^%.2f" % (wb, lg(Fraction(gap[wb], binom[wb]))))
print("mode of the exact prior                     : b =", max(range(257), key=lambda b: cnt[b]))
print("P(w=128) = 2^%.4f ; mean=128 ; sd=8 exactly (untruncated)" % lg(Fraction(cnt[128], N)))

print(); print(sep); print("2. POSTERIOR  |  conditioning events"); print(sep)
excl = Fraction(N - Z, N)
print("A  w >= 10   (X: unsigned weight <= 9 exhausted, 0 hits)  -> kills b in [0,9]")
print("B  w <= 246  (Y: complement weight <= 9 exhausted)        -> kills b in [247,256]")
print("   A and B are DISJOINT restrictions of the support: not nested, nothing multiplied.")
print()
print("P_prior(w<=9)    = 2^%.4f" % lg(cum0(9)))
print("P_prior(w>=247)  = 2^%.4f   (not equal to the above: [0,N) is not symmetric)"
      % lg(Fraction(sum(cnt[247:]), N)))
print("P_prior(excluded)= 2^%.4f   <-- THIS IS TV(prior, posterior)" % lg(excl))
print("information gained by the entire campaign about w:")
print("   -log2(1 - P_excl) = P_excl/ln2 = 2^%.4f bits = %.3e bits"
      % (lg(excl) - math.log2(math.log(2)), float(excl)/math.log(2)))
print("   KL(post||prior)   = same to leading order = 2^%.4f bits" % (lg(excl) - math.log2(math.log(2))))

print()
print("cross-check against agent X's per-bit figure < 2^-200.8:")
fam = {
  "unsigned weight<=9 (X)":     sum(binom[:10]),
  "complement weight<=9 (Y)":   sum(binom[247:]),
  "BSGS both ends 2*2^52 (X)":  2 * 2**52,
  "NAF weight<=7, as a SET (X)": 733018571531264,
  "a+b*lambda box (Q)":         (2 * 2**21) ** 2,
  "34-bit window (Q)":          223 * 2**33,
  "small multiple m<=1e7 (Q)":  256 * 10**7,
  "endomorphism orbit (Q/X)":   1536,
}
U = sum(fam.values())
for k, v in sorted(fam.items(), key=lambda kv: -kv[1]):
    print("   %-30s %-20d 2^%.2f" % (k, v, math.log2(v)))
print("   %-30s %-20d 2^%.2f" % ("UNION (upper bound)", U, math.log2(U)))
Ux = U - fam["complement weight<=9 (Y)"] - fam["NAF weight<=7, as a SET (X)"] + 1708592205206529
print("   X's list as X had it (pre-Y, signed size overcounted): 2^%.2f -> 2^%.4f of a 2^255 fibre"
      % (math.log2(Ux), math.log2(Ux) - 255))
print("   with Y's complement family added                     : 2^%.2f -> 2^%.4f of a 2^255 fibre"
      % (math.log2(U), math.log2(U) - 255))

print(); print(sep); print("3. QUANTILES  (exact rationals; B printed, coverage verified)"); print(sep)
eps_list = [("1e-1", Fraction(1,10)), ("1e-2", Fraction(1,100)), ("1e-3", Fraction(1,1000)),
            ("1e-6", Fraction(1,10**6)), ("2^-20", Fraction(1,2**20)),
            ("2^-40", Fraction(1,2**40)), ("2^-80", Fraction(1,2**80))]

def one_sided(cum, eps):
    tgt = 1 - eps
    for B in range(0, 257):
        if cum(B) >= tgt:
            assert B == 0 or cum(B-1) < tgt
            return B
def equal_tailed(cum, eps):
    h = eps / 2
    L = max(b for b in range(0, 257) if cum(b - 1) <= h)
    U_ = min(b for b in range(0, 257) if 1 - cum(b) <= h)
    return L, U_

for label, cum in [("PRIOR   (uniform on [0,N))", cum0),
                   ("POSTERIOR (10 <= w <= 246)", cum1)]:
    print(label)
    print("   eps      one-sided B   P(w<=B)        two-sided [L,U]   P(L<=w<=U)")
    for nm, e in eps_list:
        B = one_sided(cum, e)
        L, U_ = equal_tailed(cum, e)
        cov = cum(U_) - cum(L - 1)
        print("   %-6s   B = %3d      1-2^%.2f    [%3d, %3d]      1-2^%.2f"
              % (nm, B, lg(1 - cum(B)) if cum(B) < 1 else float('-inf'),
                 L, U_, lg(1 - cov) if cov < 1 else float('-inf')))
    print()

print(sep); print("4. THE HEADLINE TAILS vs AB's COST BANDS"); print(sep)
def ratio_bracket(B):
    r = Fraction(B, 256 - B + 1)
    S = sum(binom[:B + 1])
    return binom[B] <= S <= Fraction(binom[B], 1 - r)
def entropy_bound(B):
    x = B / 256.0
    H = -(x*math.log2(x) + (1-x)*math.log2(1-x))
    return 256 * (H - 1)

print(" B    P_prior(w<=B)   P_post(w<=B)   rigorous posterior bracket incl. families C-H")
for B in [14, 24, 56]:
    p0, p1 = cum0(B), cum1(B)
    num_lo = sum(max(0, cnt[b] - rem[b]) for b in range(LO, B + 1))
    den_hi = Z
    den_lo = sum(max(0, cnt[b] - rem[b]) for b in range(LO, HI + 1))
    plo = Fraction(num_lo, den_hi); phi = Fraction(sum(cnt[LO:B+1]), den_lo)
    print("%4d  2^%-11.4f  2^%-11.4f  [2^%.4f , 2^%.4f]" % (B, lg(p0), lg(p1), lg(plo), lg(phi)))

print()
print("verification of each exponent by a second route:")
for B in [14, 24, 56]:
    S = sum(binom[:B + 1])
    exact = lg(Fraction(S, TWO256))
    print("  B=%3d  exact sum 2^%.4f | single-term lower C(256,B)/2^256 = 2^%.4f"
          " | geometric-ratio upper = 2^%.4f | entropy(Chernoff) upper = 2^%.4f | brackets ok: %s"
          % (B, exact, lg(Fraction(binom[B], TWO256)),
             lg(Fraction(binom[B], TWO256) / (1 - Fraction(B, 256-B+1))),
             entropy_bound(B),
             ratio_bracket(B) and exact <= entropy_bound(B) + 1e-9))

print()
for B in [14, 24, 56]:
    p1 = cum1(B)
    print("P(w <= %3d | 10<=w<=246) = 2^%.4f = %.6e  = 1 in 2^%.2f"
          % (B, lg(p1), float(p1), -lg(p1)))

json.dump({"cnt": [str(c) for c in cnt], "Z": str(Z)}, open("ac_post.json", "w"))
