#!/usr/bin/env python3
"""Agent AG -- RED TEAM on Theorem B.  PART 2: the three attacks.

ATTACK 1  covering CODE instead of covering: does amortisation across balls buy anything?
ATTACK 2  the quantifier: is every search-based upper bound a covering of {wt>B}?
ATTACK 3  memory-aware break-even -- AB corrected the reach table for memory and left the
          break-even/crossover at their unbounded-memory values.
"""
import math
from math import comb, log2
from ag_recompute import (V128, V256, above, rep_ab, rep_tight, raw, suffix_min,
                          SUF_AB, SUF_TIGHT, ballcost, cover, rect_floor, L2, RHO)

C_VOW = 2.5

# ------------------------------------------------------------------ memory-aware primitives
def halflist(W):
    """entries in one half-list of a radius-W ball MITM (one partition)"""
    return V128[min((W + 1)//2, 128)]

def ball_time(W, logM, repf=rep_ab):
    """log2 time to search ONE radius-W ball with logM entries of memory.
       classical MITM if the half-list fits; van Oorschot-Wiener golden collision if not."""
    if W <= 0: return 0.0
    Lh = halflist(W); r = repf(W)
    if L2(Lh) <= logM:
        return L2(r*Lh)
    return L2(r*C_VOW) + 1.5*L2(Lh) - 0.5*logM

def ball_time_suffix(logM, repf=rep_ab):
    """monotone (suffix-minimised) ball time, same repair AB applied to the time-only model"""
    S = [0.0]*258
    S[257] = float('inf')
    S[256] = ball_time(256, logM, repf)
    for W in range(255, -1, -1):
        S[W] = min(ball_time(W, logM, repf), S[W + 1])
    return S

def cover_mem(B, logM, repf=rep_ab):
    """log2 cost of proving w<=B by ball covering with logM memory; capped by rho, which
       solves outright with O(1) memory and RETURNS w EXACTLY."""
    R = above(B)
    if R == 0: return (float('-inf'), None)
    S = ball_time_suffix(logM, repf)
    best = (float('inf'), None)
    for W in range(0, 257):
        t = L2(max(1.0, R/V256[W])) + S[W]
        if t < best[0]: best = (t, W)
    return best

def breakeven(logM, repf=rep_ab):
    for B in range(255, -1, -1):
        if cover_mem(B, logM, repf)[0] >= RHO:
            return B
    return -1

def crossover(logM, repf=rep_ab):
    """largest w whose low-weight ball search still costs <= rho"""
    S = ball_time_suffix(logM, repf)
    c = [w for w in range(1, 257) if S[w] <= RHO]
    return max(c) if c else 0

def reach(logT, logM, repf=rep_ab):
    S = ball_time_suffix(logM, repf)
    c = [w for w in range(1, 257) if S[w] <= logT]
    return max(c) if c else 0


if __name__ == "__main__":
    print("="*94)
    print("AG PART 2 -- ATTACKS")
    print("="*94)

    # ================================================================= ATTACK 1
    print()
    print("### ATTACK 1: covering CODE vs covering.  Can amortisation across balls help? ###")
    print()
    print("MY MODEL (different object from AB's).  Split the 256 exponent positions in half.")
    print("Because the exponents are DISTINCT POWERS OF TWO, k = k_H*2^128 + k_L holds over Z with")
    print("NO CARRIES, so wt(k) = wt(k_H) + wt(k_L) exactly.  Every meet-in-the-middle -- ball,")
    print("code, design, or ad hoc -- stores a set A of high halves and streams a set B of low")
    print("halves, and certifies precisely the COMBINATORIAL RECTANGLE A x B.  Two balls that")
    print("'share most of their half-lists' are, in this model, two rectangles with overlapping")
    print("sides -- and if they share the SAME A they merge into one rectangle A x (B1 u B2) at no")
    print("extra cost.  So the model already contains every amortisation the coordinator asked about.")
    print()
    print("FLOOR.  A rectangle of area z costs >= sqrt(z) to search by ANY algorithm (classical")
    print("MITM |A|+|B| >= 2 sqrt(z); memory-free vOW/rho on the same rectangle >= sqrt(z)).  To")
    print("prove w<=B the rectangles must cover {(a,b): wt(a)+wt(b)>B}, of size Z=|{wt>B}|.  Then")
    print("      total >= sum_i sqrt(z_i) >= sqrt(sum_i z_i) >= sqrt(Z).")
    print("Subadditivity of sqrt: the bound is INDIFFERENT to how the covering is structured, so")
    print("NO covering code, design, overlap or shared half-list can beat sqrt(Z).")
    print()
    print("SECOND, INDEPENDENT DERIVATION (generic-group counting, not combinatorics).  After m")
    print("generic queries the algorithm holds m elements sigma(a_i+b_i k); it can rule out a")
    print("scalar x only if x is a root of one of the <= C(m,2) affine collision equations or of")
    print("one of the m direct tests.  Ruling out Z scalars therefore needs m >= sqrt(2Z).  Same")
    print("floor, from a model with no rectangles in it at all.")
    print()
    print("  B     COVER (AB model)   sqrt(Z) floor    HEADROOM for any amortisation")
    for B in (250,245,240,230,220,200,180,160,152,148,140,128):
        t,_ = cover(B)
        print("  %4d    2^%8.3f       2^%8.3f        2^%.2f"%(B, t, rect_floor(B), t-rect_floor(B)))
    mx = max(cover(B)[0]-rect_floor(B) for B in range(120, 252))
    print()
    print("  MAX HEADROOM over B in [120,251]:  2^%.2f  (a factor of %.1f)"%(mx, 2**mx))
    print("  VERDICT ATTACK 1: FAILS.  AB's per-ball, no-amortisation pricing sits at most 2^%.1f"%mx)
    print("  above an amortisation-proof floor.  A covering CODE cannot move the break-even by more")
    print("  than ~3 in B.  AB's model is CONSERVATIVE, not wrong, and the conservatism is bounded.")

    # ================================================================= ATTACK 2
    print()
    print("### ATTACK 2: the quantifier.  'Every search-based upper bound is a covering of {wt>B}' ###")
    print()
    print("FALSE AS STATED, and the counterexample is the algorithm the campaign is already running.")
    print("To certify w<=B you may search {wt<=B} instead of covering {wt>B}.  A HIT hands you k,")
    print("hence a certified (indeed exact) w<=B.  That is search-based, unconditional on success,")
    print("and it is a covering of the COMPLEMENT of the set Theorem B says must be covered.")
    print()
    print("   B    Theorem-B covering cost   search-{wt<=B} cost   AB OVERPRICES BY   P(w<=B | null)")
    from math import erf, sqrt
    def pnull_le(B):
        return sum(comb(256,j) for j in range(0,B+1))/2**256
    for B in (10,14,20,30,40,52,56,64,80,100,106,110,120,128):
        cv,_ = cover(B); sm = L2(ballcost(B))
        print("  %4d       2^%8.3f            2^%8.3f          2^%7.3f        %.3e"
              %(B, cv, sm, cv-sm, pnull_le(B)))
    print()
    print("  So Theorem B's cost claim is WRONG BY UP TO 2^98 for small B.  What rescues the")
    print("  QUALITATIVE conclusion is not cost -- it is SUCCESS PROBABILITY: this certifier")
    print("  returns an upper bound only in the branch where the bound is true, and in that branch")
    print("  it has produced k itself.  Theorem B needs the hypothesis 'zero-error / works for")
    print("  every w', which is nowhere in its statement.")
    print()
    print("  THE UNIFIED LAW (mine, subsuming AB's Theorem B and matching AB's Theorem D at every B,")
    print("  not just B>=128):  cost of deciding [w<=B]  =  Theta( sqrt( min(|{w<=B}|,|{w>B}|) ) ).")
    print("     B     sqrt(min(.,.))   best known algorithm      which side")
    for B in (10,20,40,60,80,100,106,120,128,140,148,180,220,245):
        lo = V256[min(B,256)]; hi = above(B)
        floor = 0.5*L2(min(lo,hi))
        cv,_ = cover(B); sm = L2(ballcost(B))
        best = min(cv, sm, RHO)
        which = "search {w<=B}" if sm <= cv else "cover {w>B}"
        if RHO < min(cv,sm): which = "rho (solve)"
        print("  %4d      2^%8.3f       2^%8.3f              %s"%(B, floor, best, which))
    print()
    print("  AB's round-2 line 'the gap to the generic lower bound is now <= 2^3 EVERYWHERE,")
    print("  nothing large is left on the table' is FALSE outside the range AB tabulated (B>=128).")
    print("  At B=20 AB's model says 2^128.0 and the generic floor is 2^49.0 -- and the floor is")
    print("  the truth, because the low-weight sweep achieves 2^%.1f."%L2(ballcost(20)))
    print()
    print("  OTHER CANDIDATES TESTED AGAINST THE QUANTIFIER:")
    print("   (a) interval / kangaroo: decide k0 < 2^t at cost ~2^(t/2+1); a HIT proves w<=t.")
    print("       At t=128 that is 2^65 for a certified w<=128 -- 2^63 BELOW AB's curve.  It is")
    print("       search-based and is NOT a Hamming-ball covering (it is an interval, not a ball).")
    for t in (64,100,128,160,200):
        print("         t=%3d  kangaroo 2^%.1f  proves w<=%3d on a hit;  P(hit|null)=2^%.1f"
              %(t, t/2+1, t, t-256))
    print("       Excluded by success probability, NOT by Theorem B.  Same trichotomy branch as (2).")
    print("   (b) algebraic certificate (Nullstellensatz/Positivstellensatz refutation of")
    print("       {ladder eqns} ^ {sum s_i > B}).  NOT a covering, and NOT generic (it works in the")
    print("       coordinate ring of E over F_p, which Theorem D explicitly excludes).  Note the")
    print("       certificate SIZE is not the barrier at all: k0 itself is a 256-bit certificate")
    print("       verifiable by one scalar multiplication, so a proof of w<=B is always SHORT.")
    print("       The barrier can only be the cost of FINDING one, which is exactly the d_reg")
    print("       question of AB's S4 -- measured only at n=2,3.  NEITHER THEOREM COVERS THIS.")
    print("   (c) a decision procedure returning w without locating S: covered generically by")
    print("       Theorem D, NOT covered non-generically.")
    print("   (d) non-search arithmetic bounds (k0<N + digit-DP -> w<=255): already outside the")
    print("       quantifier by AB's own S7, which proves the quantifier is doing real work.")

    # ================================================================= ATTACK 3
    print()
    print("### ATTACK 3: memory-aware break-even.  AB fixed the reach table and not the headline ###")
    print()
    print("AB's S3 corrected the (time,memory)->w table with vOW (T = rep*c*L^1.5/sqrt(M)) and")
    print("published crossover w~52-64.  But 'crossover w=106' and 'break-even B=148' in the")
    print("AUTHORITATIVE SUMMARY are the UNBOUNDED-MEMORY numbers.  Re-derived with the same vOW")
    print("model AB itself adopted:")
    print()
    print("   memory      crossover w    break-even B    band [break-even, crossover] width")
    for logM in (30, 35, 40, 50, 60, 80, 10**6):
        lab = "unbounded" if logM > 1000 else "2^%d"%logM
        xo = crossover(logM); be = breakeven(logM)
        print("   %-10s      %4d           %4d              %d"%(lab, xo, be, be-xo))
    print()
    print("   -> at this box's 2^30 the break-even is B=%d, not 148.  AB's own S5 payoff table"%breakeven(30))
    print("      and the campaign's planning number both use 148/106.  In ONE model (vOW, memory")
    print("      M, rho at 2^%.2f) the correct pair at M=2^30 is (%d, %d)."%(RHO, crossover(30), breakeven(30)))
    print()
    print("   RE-CHECK of AB's S3 reach table (same model, my code):")
    print("   time \\ memory   2^30   2^35   2^40   2^50   unbounded    [AB published]")
    ab_pub = {40:(12,12,14,14,14),47:(14,16,16,18,18),58:(18,20,20,22,24),
              70:(22,24,24,26,32),80:(26,28,30,32,40),90:(32,32,34,36,48),
              126.5:(52,54,56,60,106)}
    for logT in (40,47,58,70,80,90,126.5):
        row = tuple(reach(logT,m) for m in (30,35,40,50,10**6))
        print("   2^%-6s      %4d   %4d   %4d   %4d   %6d       %s  %s"
              %(logT,row[0],row[1],row[2],row[3],row[4], ab_pub[logT],
                "MATCH" if row==ab_pub[logT] else "<-- DIFFERS"))
    print()
    print("   'w<=14 at 2^47 time / 2^30 memory':  my recomputation gives w <= %d"%reach(47,30))
    print("   'w<=52 at the rho crossover, 2^30 memory':  my recomputation gives w <= %d"%crossover(30))
