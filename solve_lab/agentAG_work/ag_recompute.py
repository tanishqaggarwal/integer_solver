#!/usr/bin/env python3
"""Agent AG -- RED TEAM on Theorem B.  PART 1: arithmetic audit, recomputed from scratch.

Nothing here is copied from ab_mem.py / ab_costfix.py.  Where I reach the same number I say so
and say whether the derivation was independent or a re-check of AB's formula.

Independent objects built here:
  * V256[j], V128[j]           cumulative Hamming volumes (exact ints)
  * rep(W)                     re-derived from first principles as 1/P[hypergeometric balance]
  * ballcost(W)                per-ball MITM cost, suffix-minimised
  * cover(B)                   min over W of (#balls) * ballcost(W)
  * RECT(B)                    MY MODEL: rectangle-covering lower bound  = sqrt(|{wt>B}|)
  * small(B)                   MY MODEL: cost of certifying w<=B by searching {wt<=B} instead
"""
import math
from math import comb, log2

def L2(x):
    return log2(x) if x > 0 else float('-inf')

# ---------------------------------------------------------------- exact volumes
V256 = [0]*257
s = 0
for j in range(257):
    s += comb(256, j); V256[j] = s
assert V256[256] == 1 << 256

V128 = [0]*129
s = 0
for j in range(129):
    s += comb(128, j); V128[j] = s
assert V128[128] == 1 << 128

def above(B):
    """|{x in F_2^256 : wt(x) > B}|"""
    return (1 << 256) - V256[min(B, 256)] if B >= 0 else 1 << 256

TOT = comb(256, 128)

# ---------------------------------------------------------------- rep(W), rederived
# A splitting system: we need, for every flip-set D with |D| = W, at least one partition of the
# 256 positions into halves L,R with |D n L| <= ceil(W/2) and |D n R| <= ceil(W/2).  A uniformly
# random 128/128 partition achieves that for a FIXED D with probability
#     P = sum over the balanced splits of the hypergeometric pmf.
# For even W the only balanced split is (W/2, W/2); for odd W both (a, W-a) with a = (W+-1)/2
# satisfy the constraint, so TWO hypergeometric terms are admissible.  AB's rep() uses only the
# a = ceil(W/2) term.  I compute BOTH and report the difference -- this is a place AB could be
# conservative (charging up to 2x too much) and I want it named.
def rep_ab(W):
    if W <= 0: return 1.0
    a = (W + 1)//2
    num = comb(W, a)*comb(256 - a*0 - W, 128 - a)   # C(256-W, 128-a)
    return TOT/num if num > 0 else float('inf')

def rep_tight(W):
    """1 / P[a random 128/128 split leaves <= ceil(W/2) on each side]  -- the honest quantity."""
    if W <= 0: return 1.0
    c = (W + 1)//2
    num = 0
    for a in range(W - c, c + 1):           # both sides <= c
        if 0 <= 128 - a <= 256 - W:
            num += comb(W, a)*comb(256 - W, 128 - a)
    return TOT/num if num > 0 else float('inf')

def raw(W, repf):
    if W <= 0: return 1.0
    return repf(W)*V128[min((W + 1)//2, 128)]

def suffix_min(repf):
    S = [0.0]*258
    S[257] = float('inf')
    S[256] = raw(256, repf)
    for W in range(255, -1, -1):
        S[W] = min(raw(W, repf), S[W + 1])
    return S

SUF_AB    = suffix_min(rep_ab)
SUF_TIGHT = suffix_min(rep_tight)

def ballcost(W, S=SUF_AB):  return S[max(0, min(W, 256))]

def cover(B, S=SUF_AB):
    """min over W of  ceil(|{wt>B}|/Vol256(W)) * ballcost(W)   -- log2."""
    R = above(B)
    if R == 0: return float('-inf')
    best = (float('inf'), None)
    for W in range(0, 257):
        nb = max(1.0, R/V256[W])
        t = L2(nb) + L2(ballcost(W, S))
        if t < best[0]: best = (t, W)
    return best

# ---------------------------------------------------------------- MY MODEL A: rectangles
# Any MITM certifies exactly a combinatorial RECTANGLE.  Split the 256 exponent positions into
# two halves; k decomposes as k = k_H*2^128 + k_L with NO carries, so wt(k) = wt(k_H)+wt(k_L).
# An MITM stores A subset F_2^128 (high halves) and streams B subset F_2^128 (low halves); it
# certifies "no solution with (k_H,k_L) in A x B".  Cost |A|+|B|, and >= sqrt(|A|*|B|) even for
# a memory-free (rho/vOW) search of the same rectangle.  To prove w<=B you must cover
# {(a,b): wt(a)+wt(b) > B}, of size |{wt>B}|.  Then
#     total cost >= sum_i sqrt(area_i) >= sqrt(sum_i area_i) >= sqrt(|{wt>B}|)
# regardless of centres, regardless of how much the balls overlap, and regardless of whether the
# covering is a code, a design or a random pile.  THIS IS THE AMORTISATION-PROOF FLOOR.
def rect_floor(B):
    return 0.5*L2(above(B))

# ---------------------------------------------------------------- MY MODEL B: search the OTHER side
# Theorem B assumes a proof of w<=B must cover {wt>B}.  It need not: search {wt<=B} instead.
# If the search HITS you have k, hence a certified w<=B (and the exact w).  If it MISSES you have
# proved w>B -- a lower bound, not an upper bound.  So this is a one-sided certifier, and it is
# not a covering of {wt>B}.  Its cost:
def small(B):
    return L2(ballcost(min(B, 256)))

RHO = 128 + math.log2(math.sqrt(math.pi/4)/math.sqrt(6))    # sqrt(pi N /4)/sqrt(6)

if __name__ == "__main__":
    print("="*94)
    print("AG PART 1 -- ARITHMETIC AUDIT OF THEOREM B, recomputed from scratch")
    print("="*94)
    print()
    print("-- A1. rho reference point (recomputed, not quoted) --")
    print("   sqrt(pi*N/4)/sqrt(6) = 2^%.4f    AB/campaign quote 2^126.5"%RHO)
    print()
    print("-- A2. rep(W): AB's single-term form vs the honest two-term form --")
    print("   W    rep_AB      rep_tight    ratio      sqrt(pi W/2)")
    for W in (2,3,9,10,11,20,21,55,106,107,127,128,255,256):
        print("  %4d  %10.4f  %10.4f  %7.4f    %8.4f"
              %(W, rep_ab(W), rep_tight(W), rep_ab(W)/rep_tight(W), math.sqrt(math.pi*W/2)))
    print("   -> for ODD W, AB's rep is up to 2.00x too LARGE (it counts only one of the two")
    print("      admissible splits).  That is CONSERVATIVE (it overcharges the attacker), so it")
    print("      cannot have produced a spurious barrier; but it is an error of up to 1 bit.")
    print()
    print("-- A3. the W=256 self-certificate --")
    for nm, S in (("AB rep", SUF_AB), ("tight rep", SUF_TIGHT)):
        print("   %-10s ballcost(256) = 2^%.6f   rep(256)=%.6f   exact 2^128? %s"
              %(nm, L2(ballcost(256, S)),
                (rep_ab if nm=="AB rep" else rep_tight)(256),
                abs(L2(ballcost(256,S))-128.0) < 1e-9))
    print("   NOTE this certificate CANNOT FAIL for a model of the form rep*Vol128(ceil(W/2))")
    print("   whose rep(256)=1, because Vol128(128)=2^128 identically.  It tests the boundary of")
    print("   rep and nothing else.  It is a weak certificate -- see AG PART 3.")
    print()
    print("-- A4. monotonicity --")
    for nm, S in (("AB rep", SUF_AB), ("tight rep", SUF_TIGHT)):
        mono = all(S[W] <= S[W+1]*(1+1e-12) for W in range(1, 256))
        rawmono = all(raw(W, rep_ab if nm=="AB rep" else rep_tight)
                      <= raw(W+1, rep_ab if nm=="AB rep" else rep_tight)*(1+1e-12)
                      for W in range(1,256))
        print("   %-10s suffix-min monotone: %s     raw (pre-suffix) monotone: %s"%(nm, mono, rawmono))
    print()
    print("-- A5. crossover and break-even, recomputed --")
    for nm, S in (("AB rep", SUF_AB), ("tight rep", SUF_TIGHT)):
        xover = max(w for w in range(1, 257) if L2(ballcost(w, S)) <= RHO)
        be = None
        for B in range(255, -1, -1):
            if cover(B, S)[0] >= RHO: be = B; break
        print("   %-10s  rho crossover w = %3d      break-even B = %3d"%(nm, xover, be))
    print("   AB published: crossover 106, break-even 148.  Z: 106 / 149.")
    print()
    print("-- A6. the covering-cost curve, and where it saturates --")
    print("     B    |{wt>B}|   bestW   #balls    cost/ball    COVER      sqrt(|{wt>B}|)  gap")
    for B in (250,245,240,230,220,200,180,160,152,148,140,128,106,80,40,20,0):
        t, W = cover(B)
        R = above(B); nb = max(1.0, R/V256[W])
        print("  %4d   2^%6.1f   %4d   2^%6.1f   2^%7.2f   2^%7.3f   2^%7.3f     2^%.2f"
              %(B, L2(R), W, L2(nb), L2(ballcost(W)), t, rect_floor(B), t-rect_floor(B)))
    print()
    print("   *** The covering curve is CAPPED at 2^128.000 for every B <= 148: the minimising W")
    print("   *** is 256 (one ball = the whole space).  So the whole family of statements")
    print("   *** 'prove w<=B' for B in [0,148] costs between 2^%.1f and 2^128.000 -- a band of"%RHO)
    print("   *** 2^%.2f.  There is no cliff at 148; the curve is FLAT below it."%(128.0-RHO))
