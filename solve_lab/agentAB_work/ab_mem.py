#!/usr/bin/env python3
"""Agent AB round 3.  Fixes Z's four items, then the memory-aware MITM costing.

FIXED ball cost.  Two defects in round 2's  sqrt(W)*Vol_128(W//2):
  (i)  floor vs ceil: an odd-weight flip set splits (ceil(W/2), floor(W/2)); the cost is set by
       the LARGER half, so it must be Vol_128(ceil(W/2)).  My scans stepped even W only, which
       hid it; it fires on first reuse.
  (ii) the sqrt(W) partition factor is WRONG AT THE BOUNDARY.  The number of random 128/128
       partitions needed is the reciprocal of the probability that one cuts a fixed W-set evenly:
            P(W) = C(W,W/2)*C(256-W,128-W/2) / C(256,128)     (exact, hypergeometric)
       P -> sqrt(2/(pi W)) for small W (recovering sqrt(W)), but P(256) = 1 exactly, so at W=256
       the factor is 1, not 16.  That is why round 2's certificate returned 2^132.0 where the
       true answer is 2^128.0 -- and I printed the miss and moved on.  Same failure mode I had
       just retracted.  Fixed here and re-certified as an EXACT match.
"""
import math
from math import comb, log2
def L2(x): return log2(x) if x>0 else float('-inf')
N=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
V128=[0]*129; s=0
for j in range(129): s+=comb(128,j); V128[j]=s
TOT=comb(256,128)
def rep(W):
    """expected number of random half/half partitions until one splits a fixed W-set evenly"""
    if W<=0: return 1.0
    a=(W+1)//2
    num=comb(W,a)*comb(256-W,128-a)
    return TOT/num if num>0 else float('inf')
def raw_cost(W):
    if W<=0: return 1.0
    return rep(W)*V128[min((W+1)//2,128)]
# A radius-W ball can always be searched by the radius-W' procedure for any W' >= W, so the true
# cost is the SUFFIX MINIMUM.  This also repairs the only non-monotonicity, at the very top where
# rep(255)=2 makes the raw cost of radius 255 exceed that of radius 256.
_SUF=[0.0]*258; _SUF[257]=float('inf'); _SUF[256]=raw_cost(256)
for _W in range(255,-1,-1): _SUF[_W]=min(raw_cost(_W),_SUF[_W+1])
def ball_cost(W): return _SUF[max(0,min(W,256))]
def Vol(W): return sum(comb(256,j) for j in range(W+1))
def region(B): return sum(comb(256,j) for j in range(B+1,257))

print("== Z item 2: the W=256 certificate, re-run with the fixed cost ==")
print("   required (one ball = whole key space, birthday bound) : 2^128.0")
print("   round-2 model  sqrt(W)*Vol_128(W//2)                  : 2^%.4f   <- MISSED BY 2^4, I printed it"%L2(16*V128[128]))
print("   fixed model    rep(W)*Vol_128(ceil(W/2))              : 2^%.4f   rep(256)=%.4f"%(L2(ball_cost(256)),rep(256)))
print("   EXACT match:", abs(L2(ball_cost(256))-128.0)<1e-9)
print()
print("== Z item 1+3: monotonicity, and the floor/ceil bug ==")
mono=all(ball_cost(W)<=ball_cost(W+1)*1.0000001 for W in range(1,256))
print("   fixed ball_cost (suffix-min) monotone non-decreasing on W=1..256 :", mono)
print("   round-1 C(256,W/2) monotone :", all(comb(256,W//2)<=comb(256,(W+1)//2) for W in range(1,256)),
      "  <- Z's third disqualifier, which I had not stated")
print("   odd radii, floor vs ceil (the latent bug):")
for W in (9,15,35,55,75,103,127):
    print("     W=%3d  floor 2^%6.2f   ceil(correct) 2^%6.2f   underpriced by %.1f bits"
          %(W,L2(rep(W)*V128[W//2]),L2(raw_cost(W)),L2(raw_cost(W))-L2(rep(W)*V128[W//2])))
print("   rep(W) sanity: rep(10)=%.2f vs sqrt(pi*10/2)=%.2f ; rep(256)=%.4f"%(rep(10),math.sqrt(math.pi*10/2),rep(256)))

print()
print("== Z item 4: THE HEADLINE, RESTATED IN ONE MODEL ==")
print(" Corrected constant (Z is right, my m^2/(2 min) was a factor 2 too tight):")
print("   Adv <= C(m,2)*(1/|D0| + 1/|D1|) <= (m^2/2)*(2/min) = m^2/min   ->  m >= sqrt(eps*min)")
mn=0.0
# exact digit-DP count of {k < N : popcount(k) <= B}, so the universe is [0,N) not [0,2^256)
def count_le(B):
    bits=[(N>>i)&1 for i in range(256)]
    from functools import lru_cache
    @lru_cache(None)
    def f(i,used,tight):
        if used>B: return 0
        if i<0: return 1
        tot=0
        for b in (0,1):
            if tight and b>bits[i]: continue
            tot+=f(i-1,used+b, tight and b==bits[i])
        return tot
    r=f(255,0,True); f.cache_clear(); return r
assert count_le(256)==N+1, "digit-DP self-check failed"   # DP counts [0,N] inclusive
print("   digit-DP self-check  #{0<=k<=N} == N+1 :", count_le(256)==N+1)
print("   popcount(N)=%d > 128, so N itself is not in D0 and the [0,N) counts are exact"%bin(N).count('1'))
a=count_le(128); b=N-a; mn=min(a,b)
print("   B=128 over the true universe [0,N): |D0|=2^%.2f |D1|=2^%.2f min=2^%.2f"%(L2(a),L2(b),L2(mn)))
base=0.5*L2(mn)-0.5                       # m >= sqrt(eps*min), eps=1/2
print("   generic bound, NO automorphisms, eps=1/2 :  m >= 2^%.2f"%base)
aut=base-0.5*log2(6)
print("   same model as the 2^126.5 solve (encoding + GLV: order-6 orbit collapses, AUT=6")
print("   equations per pair, degrading the bound by sqrt(6)) :  m >= 2^%.2f"%aut)
print()
print("   STRUCK:  'm >= 2^127.5 vs solving 2^126.5, deciding is as hard as solving'")
print("            -- read literally that says deciding is HARDER than solving, which is")
print("               impossible: any solver decides.  The two numbers were in different models.")
print("   STANDS:  in ONE model, deciding w<=128 needs m >= 2^%.1f and solving costs 2^126.5."%aut)
print("            The bound points the RIGHT way (2^%.1f <= 2^126.5, gap 2^%.1f), and the"%(aut,126.5-aut))
print("            qualitative claim is untouched: NO generic shortcut for the weight predicate,")
print("            same cost as solving to within the same sqrt(6) the solver itself uses.")

print()
print("== Theorem B recomputed with the fixed cost, and reconciled with Z ==")
best=None
for W in range(1,257):
    if L2(ball_cost(W))<126.5: best=W
print("   largest affordable complement radius W = %d  ->  ceiling  w <= %d"%(best,255-best))
be=None
for B in range(255,100,-1):
    t=min(L2(max(1.0,region(B)/Vol(W)))+L2(ball_cost(W)) for W in range(1,257))
    if t>=126.5: be=B; break
print("   break-even B = %d (first B whose covering proof costs >= 2^126.5)"%be)
wx=max(w for w in range(1,257) if L2(ball_cost(w))<=126.5)
print("   rho crossover: largest w whose ball search costs <= 2^126.5  =  %d"%wx)
print("   RECONCILIATION with Z (crossover 106, break-even 149) and with my RESUME (104/107/148):")
print("     - my RESUME's 104 and my script's 107 were DIFFERENT quantities (budget-table reach at")
print("       2^126 vs largest affordable radius); that was my sloppiness, not a disagreement.")
print("     - with the boundary + ceil fixes I now get crossover %d, break-even %d."%(wx,be))
print("     - residual gap to Z is from the partition factor; Z's model and mine differ there.")
print("       Both agree the qualitative claim is unchanged.  Quoting a RANGE is honest: the")
print("       time-only crossover is w ~ 104-107 and the break-even ceiling B ~ 148-149.")

print()
print("="*92)
print("TASK 2 -- MEMORY-AWARE COSTING.  Achievable w as a function of (time, memory).")
print("="*92)
print(" Time-only figures assume memory = the half-list, which is the assumption that fails.")
print(" With M < L the right algorithm is NOT chunked rescanning but van Oorschot-Wiener")
print(" parallel collision search for the GOLDEN collision:")
print("     chunked passes :  T = rep * L^2 / M")
print("     vOW golden     :  T = rep * c * L^1.5 / sqrt(M),  c ~ 2.5")
print("   ratio vOW/chunked = c*sqrt(M/L) < 1 for all M < L, so vOW dominates.  Use it.")
print("   Solving for the reach:   L <= ( T^2 * M / (rep^2 * c^2) )^(1/3)   when M <= L.")
C_VOW=2.5
def reach(logT, logM):
    best=0
    for w in range(1,257):
        Lh=V128[min((w+1)//2,128)]; r=rep(w)
        if math.log2(Lh)<=logM:  t=L2(ball_cost(w))               # memory suffices
        else:                    t=L2(r*C_VOW)+1.5*L2(Lh)-0.5*logM
        if t<=logT: best=w
    return best
def halflist(w): return L2(V128[min((w+1)//2,128)])
print()
print("  reachable w   |  memory (entries):  2^30 (this box, 12GB)  2^35   2^40   2^50   unbounded")
for logT in (30,40,47,58,70,80,90,126.5):
    row=[reach(logT,m) for m in (30,35,40,50,10**6)]
    print("   time 2^%-5.1f  |                     %4d   %13d %6d %6d %10d"
          %(logT,row[0],row[1],row[2],row[3],row[4]))
print()
print("  THIS BOX (time 2^47, memory 2^30): reach w <= %d.  Z is right that my time-only 18 was"%reach(47,30))
print("  fiction -- the half-list at w=18 is 2^%.1f entries = 2^%.1f times this box's memory."
      %(halflist(18), halflist(18)-30))
print()
print("  THE CROSSOVER WITH RHO, memory-aware (rho needs O(1) memory, so this is the real question):")
for logM in (30,35,40,50,60,10**6):
    w=reach(126.5,logM)
    lab="unbounded" if logM>1000 else "2^%d"%logM
    print("    memory %-10s ->  MITM beats rho (2^126.5) only for w <= %3d"%(lab,w))
print()
print("  => the honest crossover is w ~ %d-%d for ANY realistic memory, not 106."%(reach(126.5,30),reach(126.5,60)))
print("     Memory enters as a CUBE ROOT, so the crossover is remarkably insensitive to it:")
print("     20 orders of magnitude of memory (2^30 -> 2^60) moves it by only %d."
      %(reach(126.5,60)-reach(126.5,30)))
print("     NOTE the reconciliation: the campaign's original w ~ 56 crossover was nearly right,")
print("     but for the wrong reason -- its time model was pessimistic and it ignored memory,")
print("     and the two errors happened to cancel.")
print()
print("  Disk: this box has ~30 GB free ~ 2^31 entries, and vOW's access pattern is random, so")
print("  disk-backed tables run at seek rate (~10^2/s) not memory rate (~10^8/s): a 2^20 slowdown")
print("  against a 2^1 memory gain.  Disk is NOT a way out; the reach stays w <= %d."%reach(47,30))
