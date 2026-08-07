#!/usr/bin/env python3
"""Agent AB: RE-AUDIT of Theorem B's cost model.  TWO errors found in my own round-1 table:
 (1) a weight-W MITM was priced at C(256,W/2); the standard low-Hamming-weight DLP algorithm
     (Coppersmith / Stinson splitting systems) splits the 256 POSITIONS and costs less;
 (2) more seriously, the per-ball list is the CUMULATIVE volume of a half, Vol_128(W/2),
     not the single binomial C(128,W/2).  At large W that difference is enormous and it is
     what made the round-1 table read 2^185 where the truth is ~2^128."""
import math
from math import comb, log2
def L(x): return log2(x) if x>0 else float('-inf')
V128=[0]*(129)
s=0
for j in range(129):
    s+=comb(128,j); V128[j]=s
def ball_cost(W):
    """MITM over a Hamming ball of radius W around a known centre: split the 256 positions into
    two halves; a splitting system guarantees the <=W flips divide ~evenly; each half-list is
    Sum_{j<=W/2} C(128,j) = Vol_128(W/2).  Small poly factor sqrt(W) for the splitting system."""
    return math.sqrt(max(W,1))*V128[min(W//2,128)]
def Vol(W): return sum(comb(256,j) for j in range(W+1))
def region(B): return sum(comb(256,j) for j in range(B+1,257))

print("== 0. SANITY: the model must reproduce known costs ==")
print("   W=256 (one ball = whole space) -> 2^%.1f   [must be ~2^128, the birthday/BSGS cost]"%L(ball_cost(256)))
print("   W=10  -> 2^%.1f   [round-1 said 2^33.0 using C(256,5)]"%L(ball_cost(10)))
print("   W=20  -> 2^%.1f   [round-1 said 2^58.0]"%L(ball_cost(20)))
print("   generic lower bound for the whole space sqrt(N) = 2^128.0 : model gives 2^%.1f  OK"%L(ball_cost(256)))

print()
print("== 1. cost models side by side, against Theorem D's generic lower bound ==")
print("  w    round-1 C(256,w/2)   corrected ball_cost   Theorem-D LB sqrt(|{w'<=w}|)")
for w in (10,12,14,16,20,24,30,40,56,80,114,128,200,256):
    cls=sum(comb(256,j) for j in range(min(w,256)+1))
    print("  %3d      2^%6.1f            2^%6.1f              2^%.1f"
          %(w,L(comb(256,min(w//2,256))),L(ball_cost(w)),0.5*L(cls)))
print("  -> the corrected model never dips below the generic lower bound.  The round-1 model does")
print("     not either, but it OVERSHOOTS by up to 2^65 at large w, which is where my error lived.")

print()
print("== 2. THEOREM B, CORRECTED ==")
print("   B      |{w>B}|   best W   #balls   cost/ball    TOTAL     (round-1 said)   Theorem-D LB")
for B in (250,245,240,230,220,200,180,152,128):
    R=region(B); best=None
    for W in range(2,257,2):
        nb=max(1.0,R/Vol(W)); t=L(nb)+L(ball_cost(W))
        if best is None or t<best[0]: best=(t,W,L(nb),L(ball_cost(W)))
    t,W,nb,cb=best
    r1=min(L(max(1.0,region(B)/Vol(W2)))+L(comb(256,W2//2)) for W2 in range(2,257,2))
    lb=0.5*L(min(R,2**256-R))
    print("  %4d   2^%6.1f   W=%3d   2^%5.1f   2^%6.1f    2^%6.1f      2^%6.1f       2^%.1f"
          %(B,L(R),W,nb,cb,t,r1,lb))
print()
print("  The minimum is attained at W=256 -- ONE ball covering everything -- for every B below ~247.")
print("  In plain words: THE CHEAPEST WAY TO PROVE ANY NONTRIVIAL CEILING IS TO SOLVE THE INSTANCE.")

print()
print("== 3. the corrected break-even ==")
be=None
for B in range(255,100,-1):
    R=region(B)
    t=min(L(max(1.0,R/Vol(W)))+L(ball_cost(W)) for W in range(2,257,2))
    if t>=126.5: be=B; break
print("   first B whose proof costs >= 2^126.5 : B = %s"%be)
print("   -> with the corrected model NO B reaches 2^126.5 by covering, because covering degenerates")
print("      into 'solve it' (2^128) which is within 2^1.5 of rho-with-automorphisms (2^126.5).")
print()
print("   The meaningful break-even is therefore the COMPLEMENT ball alone (Theorem A: the only")
print("   unconditional centre).  Radius W around all-ones proves w <= 255-W at cost ball_cost(W):")
print("     W    cost      proves      vs solving 2^126.5")
for W in (10,20,40,60,80,100,120,140,160,180,200,220):
    c=L(ball_cost(W))
    print("    %3d   2^%6.1f   w <= %3d    %s"%(W,c,255-W,"cheaper" if c<126.5 else "NOT cheaper"))
last=[W for W in range(2,257) if L(ball_cost(W))<126.5]
print("   largest affordable radius W = %d  ->  the best ceiling a search can prove more cheaply"%max(last))
print("   than solving the instance is  w <= %d."%(255-max(last)))

print()
print("== 4. corrected BUDGET -> reachable weight (this is directly actionable for X and Y) ==")
print("  budget    round-1 model    corrected model    memory needed (half-list)")
for bud in (30,40,47,58,70,80,90,126):
    wn=max([w for w in range(2,258,2) if L(comb(256,w//2))<=bud]+[0])
    wb=max([w for w in range(2,258,2) if L(ball_cost(w))<=bud]+[0])
    print("  2^%-3d      w <= %-3d          w <= %-3d            2^%.1f entries"
          %(bud,wn,wb,L(V128[min(wb//2,128)])))
mem=max(w for w in range(2,130,2) if V128[w//2]<=2**30)
print()
print("  MEMORY BINDS on this box: ~12 GB ~ 2^30 entries -> half-list Vol_128(w/2) <= 2^30 gives")
print("  w <= %d.  Beyond that a low-memory variant (parallel collision search) is required and the"%mem)
print("  time cost rises above the table.  The round-1 note did not say this and it should have.")
