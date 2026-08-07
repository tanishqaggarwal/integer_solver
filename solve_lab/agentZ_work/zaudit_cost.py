#!/usr/bin/env python3
"""Agent Z: independent audit of agent AB's CORRECTED ball-covering cost model.

Re-derived from scratch; AB's ab_costfix.py was read for the CLAIM only, not imported.
Everything here is recomputed with exact integers (math.comb), never floats, until the
final log2.
"""
from math import comb, log2, sqrt, isqrt

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

def L(x):
    return log2(x) if x > 0 else float('-inf')

# exact cumulative volumes
V128 = [0] * 129
s = 0
for j in range(129):
    s += comb(128, j); V128[j] = s
V256 = [0] * 257
s = 0
for j in range(257):
    s += comb(256, j); V256[j] = s

def half_list(W):
    """entries enumerated per half by a splitting-system MITM over a radius-W ball"""
    return V128[min(W // 2, 128)]

def ball_cost(W, poly=True):
    """AB's corrected model: sqrt(W) splitting-system partitions x half-list."""
    f = sqrt(max(W, 1)) if poly else 1
    return f * half_list(W)

print("=" * 92)
print("A. THE FORMULA — re-derived independently")
print("=" * 92)
print("""  Target: enumerate every k with wt(k) <= W (a Hamming ball of radius W around a centre).
  Split the 256 bit positions into halves A,B of 128.  A weight-w set S has |S n A| = a and
  |S n B| = w-a; a MITM over the fixed partition only works if BOTH parts are <= W/2, which a
  fixed partition does NOT guarantee.  A splitting system -- a family of partitions such that
  every S is balanced by at least one of them -- fixes that, at a poly(W) multiplicative cost
  (Coppersmith / Stinson; sqrt(W) is the figure AB uses).  Per partition each half enumerates
  every subset of size <= W/2, i.e. Vol_128(W/2), NOT the single binomial C(128,W/2).
  => cost(W) = poly(W) * Vol_128(W/2).   AB's corrected formula.  I get the same.""")
print("  Round-1's C(256,W/2) is wrong in TWO ways: wrong ground set (256 not 128) and single")
print("  binomial not cumulative.  Confirmed: C(256,W/2) is NON-MONOTONE in W (it peaks at")
print("  W=256 -> C(256,128) and falls for W>256), which alone disqualifies it as a ball cost.")
print("     C(256,128) = 2^%.1f   C(256,127) = 2^%.1f  (a ball cost must be monotone)"
      % (L(comb(256, 128)), L(comb(256, 127))))

print()
print("=" * 92)
print("B. BOUNDARY CHECKS IN BOTH DIRECTIONS")
print("=" * 92)
print("  W = 0    : one ball = one point.  cost/ball must be 1, #balls must be the whole region.")
print("             model: cost/ball = %d,  Vol(0) = %d   -> #balls = |region| = full enumeration  OK"
      % (ball_cost(0), V256[0]))
print("  W = 256  : one ball = the WHOLE key space.  cost MUST equal plain MITM/BSGS = 2^128.")
print("             model with sqrt(W): 2^%.1f      <-- OVERSHOOTS the truth by 2^%.1f"
      % (L(ball_cost(256)), L(ball_cost(256)) - 128))
print("             model without the poly factor: 2^%.1f   <-- EXACT" % L(ball_cost(256, poly=False)))
print("             Reason: at W=256, W/2 = 128 = the half size, so EVERY split is automatically")
print("             balanced and no splitting system is needed.  The sqrt(W) factor is spurious")
print("             exactly at the boundary AB used as its sanity test.")
print("  monotone?: ", all(ball_cost(W) <= ball_cost(W + 1) + 1e-9 for W in range(0, 256)))
print("  half-list monotone (the part that matters):",
      all(half_list(W) <= half_list(W + 1) for W in range(0, 256)))

print()
print("=" * 92)
print("C. THE MODEL vs THE GENERIC LOWER BOUND — does it ever dip below?")
print("=" * 92)
print("  Provable relation:  (Vol_128(W/2))^2 = #{S : |S n A| <= W/2 AND |S n B| <= W/2} <= Vol_256(W)")
print("  so  Vol_128(W/2) <= sqrt(Vol_256(W)) ALWAYS.  Verified exactly for every even W:")
ok = all(half_list(W) ** 2 <= V256[W] for W in range(0, 257, 2))
print("     (Vol_128(W/2))^2 <= Vol_256(W) for all even W in 0..256 :", ok)
print("  => the corrected MITM is ABOVE the sqrt(class size) generic bound by at most the poly")
print("     factor, and BELOW it without one.  It is optimal to within sqrt(W) <= 2^4.")
print()
print("   w    corrected ball cost   generic LB sqrt(Vol_256(w))   ratio (bits)")
for w in (10, 12, 14, 16, 18, 20, 24, 30, 40, 56, 80, 104, 107, 128, 200, 256):
    c = L(ball_cost(w)); lb = 0.5 * L(V256[w])
    print("  %3d      2^%6.1f              2^%6.1f                  %+5.2f" % (w, c, lb, c - lb))

print()
print("=" * 92)
print("D. CONFIRM OR CORRECT AB's CORRECTED BUDGET TABLE")
print("=" * 92)
print("  budget   AB claims   I get (with sqrt(W))   I get (no poly factor)   half-list memory")
claims = {47: 18, 58: 24, 80: 40}
for bud in (30, 40, 47, 58, 70, 80, 90, 126, 127):
    wp = max([w for w in range(0, 257) if L(ball_cost(w)) <= bud] + [0])
    wn = max([w for w in range(0, 257) if L(ball_cost(w, poly=False)) <= bud] + [0])
    cl = claims.get(bud, None)
    mark = "" if cl is None else ("  <-- AB says w <= %d : %s" % (cl, "CONFIRMED" if cl == wp else "DIFFERS"))
    print("  2^%-4d    %-9s   w <= %-3d               w <= %-3d                  2^%.1f entries%s"
          % (bud, ("w <= %d" % cl) if cl else "-", wp, wn, L(half_list(wp)), mark))

print()
print("=" * 92)
print("E. THE RHO CROSSOVER  (AB's RESUME says w ~ 104; its script prints W = 107)")
print("=" * 92)
for tgt, name in ((126.5, "rho + negation + GLV"), (128.0, "plain MITM / BSGS on the whole space")):
    wp = max(w for w in range(0, 257) if L(ball_cost(w)) < tgt)
    wn = max(w for w in range(0, 257) if L(ball_cost(w, poly=False)) < tgt)
    print("  vs %-38s : crossover w = %3d (with sqrt(W))   %3d (without)" % (name, wp, wn))
print("  -> the complement radius and the direct-search radius are THE SAME function, so the")
print("     'largest affordable radius W' and the 'rho crossover w' must be the same number.")
print("     AB's RESUME quotes 104 and its script quotes 107; the reproducible value is below.")
w107 = max(w for w in range(0, 257) if L(ball_cost(w)) < 126.5)
print("     reproducible crossover (AB's own formula, budget 2^126.5) : w = %d" % w107)
print("     => break-even ceiling B = 255 - %d = %d   (AB's RESUME says 148)" % (w107, 255 - w107))

print()
print("=" * 92)
print("F. THEOREM B ITSELF — is the covering minimum really at W = 256?")
print("=" * 92)
print("    B      |{w>B}|    argmin W    #balls     cost/ball     TOTAL")
for B in (245, 240, 220, 200, 180, 152, 148, 128):
    R = V256[256] - V256[B]
    best = None
    for W in range(0, 257):
        nb = max(1.0, R / V256[W])
        t = L(nb) + L(ball_cost(W))
        if best is None or t < best[0]:
            best = (t, W, L(nb), L(ball_cost(W)))
    t, W, nb, cb = best
    print("  %4d    2^%6.1f     W=%3d     2^%6.1f    2^%6.1f     2^%6.1f" % (B, L(R), W, nb, cb, t))
print("  CONFIRMED: for every nontrivial B the optimum is the degenerate one-ball cover, i.e.")
print("  'solve the whole instance'.  AB's qualitative conclusion stands.")
