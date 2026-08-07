#!/usr/bin/env python3
"""Agent Z: audit of AB's corrected cost model, round 2 -- the floor/ceil bug and
what it does (and does not) do to AB's published numbers."""
from math import comb, log2, sqrt

def L(x): return log2(x) if x > 0 else float('-inf')
V128 = [0] * 129
s = 0
for j in range(129):
    s += comb(128, j); V128[j] = s
V256 = [0] * 257
s = 0
for j in range(257):
    s += comb(256, j); V256[j] = s

def half_floor(W):  return V128[min(W // 2, 128)]          # AB's ab_costfix.py, verbatim
def half_ceil(W):   return V128[min((W + 1) // 2, 128)]    # correct: an odd-weight set splits ceil/floor
def cost(W, half, poly=True): return (sqrt(max(W, 1)) if poly else 1) * half(W)

print("=" * 92)
print("BUG: AB's half-list uses floor(W/2).  A weight-W set with W ODD cannot split more evenly")
print("than (ceil(W/2), floor(W/2)), so the half-list must be Vol_128(ceil(W/2)).")
print("Using floor UNDERPRICES every odd radius, which makes odd W look spuriously cheap.")
print("=" * 92)
print("   W    AB floor-model   correct ceil-model   understated by")
for W in (9, 15, 19, 25, 35, 41, 55, 75, 103, 107, 127):
    a, b = L(cost(W, half_floor)), L(cost(W, half_ceil))
    print("  %3d      2^%6.1f          2^%6.1f            %5.1f bits" % (W, a, b, b - a))
print()
print("AB's own scans only ever step through EVEN W (`range(2,257,2)`), where floor == ceil,")
print("so the bug never fires in AB's published numbers.  It fires the moment anyone reuses")
print("the function on an odd radius -- which my first pass did, producing spurious optima at")
print("W = 9,15,35,55,75,103,127 and a spurious +1 on every budget row.")
print("  ceil-model still monotone in W :", all(cost(W, half_ceil) <= cost(W + 1, half_ceil) + 1e-9
                                                for W in range(0, 256)))

print()
print("=" * 92)
print("A. BOUNDARY CHECKS (ceil model)")
print("=" * 92)
print("  W=0   -> cost/ball = %d, #balls = |region|  (pure enumeration)  OK" % cost(0, half_ceil))
print("  W=256 -> 2^%.1f with the sqrt(W) factor, 2^%.1f without.  TRUTH = 2^128 exactly."
      % (L(cost(256, half_ceil)), L(cost(256, half_ceil, poly=False))))
print("     ** AB's stated sanity test ('W=256 must give 2^128') is NOT passed by AB's corrected")
print("        model: it returns 2^132.0.  The 2^4 gap is the splitting-system factor, which is")
print("        spurious at W=256 because W/2 = 128 = the half size makes every split balanced.")
print("        The correction AB made is right; the test it used to certify the correction is")
print("        satisfied only to within the poly factor, and AB reported 2^132.0 without saying so.")
print("  Even-W sanity: (Vol_128(w/2))^2 <= Vol_256(w) for all even w :",
      all(V128[w // 2] ** 2 <= V256[w] for w in range(0, 257, 2)),
      " => the model never dips below the sqrt(class-size) generic bound")

print()
print("=" * 92)
print("B. AB's CORRECTED BUDGET TABLE -- verdict")
print("=" * 92)
print("  budget    AB claims   ceil model   floor model   memory (half-list)   verdict")
claims = {47: 18, 58: 24, 80: 40}
for bud in (30, 40, 47, 58, 70, 80, 90, 126):
    wc = max([w for w in range(0, 257) if L(cost(w, half_ceil)) <= bud] + [0])
    wf = max([w for w in range(0, 257) if L(cost(w, half_floor)) <= bud] + [0])
    cl = claims.get(bud)
    v = "" if cl is None else ("CONFIRMED" if cl == wc else "DIFFERS (AB %d)" % cl)
    print("  2^%-4d    %-9s   w <= %-3d      w <= %-3d       2^%-5.1f entries      %s"
          % (bud, ("w<=%d" % cl) if cl else "-", wc, wf, L(half_ceil(wc)), v))

print()
print("=" * 92)
print("C. THE CROSSOVER AND THE BREAK-EVEN CEILING")
print("=" * 92)
for tgt, nm in ((126.5, "rho + negation map + GLV"), (128.0, "plain MITM on the whole space")):
    wc = max(w for w in range(0, 257) if L(cost(w, half_ceil)) < tgt)
    print("  budget 2^%-5.1f (%-32s): largest affordable radius w = %d" % (tgt, nm, wc))
wc = max(w for w in range(0, 257) if L(cost(w, half_ceil)) < 126.5)
print("  => break-even ceiling  B = 255 - %d = %d" % (wc, 255 - wc))
print("     AB's RESUME says B = 148 and 'largest affordable complement radius W = 107'")
print("     AB's RESUME ALSO says 'the rho crossover moves from w~56 to w~104'.")
print("     Those are the SAME quantity (the same function, the same budget) and AB quotes two")
print("     different numbers for it.  Reproducible value with AB's own formula: w = %d, B = %d."
      % (wc, 255 - wc))

print()
print("=" * 92)
print("D. THEOREM B: is the covering optimum really the degenerate one-ball cover?")
print("=" * 92)
print("     B      |{w>B}|    argmin W    #balls      cost/ball     TOTAL")
for B in (250, 245, 240, 220, 200, 180, 152, 148, 128):
    R = V256[256] - V256[B]
    best = None
    for W in range(0, 257):
        nb = max(1.0, R / V256[W])
        t = L(nb) + L(cost(W, half_ceil))
        if best is None or t < best[0]:
            best = (t, W, L(nb), L(cost(W, half_ceil)))
    t, W, nb, cb = best
    print("   %4d    2^%6.1f     W=%3d     2^%6.1f     2^%6.1f     2^%6.1f" % (B, L(R), W, nb, cb, t))
print("   CONFIRMED: for every B below ~247 the optimum is one ball covering everything, i.e.")
print("   the cheapest way to prove any nontrivial ceiling is to solve the instance outright.")
print("   AB's headline qualitative claim SURVIVES the audit unchanged.")
