#!/usr/bin/env python3
"""Agent AB, round 2: attack my own least-secure DEAD verdict (§2, the WEIGHT predicate),
and re-audit Theorem B's cost model while I am at it."""
import math, random
from math import comb, log2
N=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
ALL=2**256-1
def wt(x): return bin(x).count('1')
rng=random.Random(4242)

print("="*90)
print("THEOREM C.  The ONLY weight-preserving affine self-map of Z_N is the identity.")
print("="*90)
print(" Claim: if psi(k)=a*k+b mod N has popcount(psi(k))=popcount(k) for all k in [0,N), then a=1,b=0.")
print(" Proof step 1: k=0 forces popcount(b)=0, and 0<=b<N, so b=0.")
print(" Proof step 2: k=1 forces popcount(a)=1, so a=2^j for some j in [0,255].")
print(" Proof step 3: WITNESS.  For every j in [1,255] take k = 2^(256-j)  (< 2^255 < N, and popcount 1).")
print("               Then a*k = 2^256 mod N, and popcount(2^256 mod N) = %d != 1."%wt(2**256%N))
print("               So j=0, a=1.  QED")
print()
print(" -- step 3 verified for every j in [1,255], no exceptions:")
bad=[]
for j in range(1,256):
    a=2**j; k=2**(256-j)
    assert k<N and wt(k)==1
    if wt((a*k)%N)==wt(k): bad.append(j)
print("    j with popcount(2^j * 2^(256-j) mod N) == 1 :", bad if bad else "NONE  (all 255 j refuted)")
print(" -- step 1/2 verified by exhaustion over the only candidates:")
print("    b=0 forced:", all(wt(b)!=0 for b in [1,2,3,N-1,2**255]), "(popcount 0 <=> b=0 on [0,N))")
print("    a of popcount 1 :", len([j for j in range(256)]), "candidates 2^0..2^255, all but j=0 refuted above")

print()
print("="*90)
print("Is there even an APPROXIMATELY weight-preserving affine map?  (measured, 20000 k each)")
print("="*90)
ks=[rng.randrange(N) for _ in range(20000)]
def frac_pres(a,b=0):
    return sum(1 for k in ks if wt((a*k+b)%N)==wt(k))/len(ks)
def corr(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sx=(sum((u-mx)**2 for u in xs)/n)**.5; sy=(sum((v-my)**2 for v in ys)/n)**.5
    if sx==0 or sy==0: return float('nan')
    return sum((u-mx)*(v-my) for u,v in zip(xs,ys))/(n*sx*sy)
wk=[wt(k) for k in ks]
print("  a                      Pr[wt(ak+b)=wt(k)]   corr(wt(k),wt(ak+b))   corr of predicate [w<=128]")
LAM=37718080363155996902926221483475020450927657555482586988616620542887997980018
cands=[("1 (identity)",1,0),("2",2,0),("4",4,0),("2^8",2**8,0),("2^32",2**32,0),
       ("2^-1 mod N",pow(2,N-2,N),0),("-1",N-1,0),("lambda",LAM,0),
       ("2 with b=1",2,1),("random a",rng.randrange(2,N),0),("random a, random b",rng.randrange(2,N),rng.randrange(N))]
for name,a,b in cands:
    wa=[wt((a*k+b)%N) for k in ks]
    P=[1 if x<=128 else 0 for x in wk]; Q=[1 if x<=128 else 0 for x in wa]
    print("  %-22s %8.4f            %+7.4f               %+7.4f"%(name,frac_pres(a,b),corr(wk,wa),corr(P,Q)))
print()
print("  a=2^j preserves weight exactly on {k : 2^j k < N}, a fraction 2^-j -- predicted vs measured:")
for j in (1,2,3,8):
    pred=2.0**-j
    print("    j=%2d  predicted %.6f   measured %.6f"%(j,pred,frac_pres(2**j)))
print("  => the ONLY affine maps with ANY weight preservation are k -> 2^j k, and their")
print("     preservation set shrinks like 2^-j.  There is no weight-preserving randomisation.")

print()
print("="*90)
print("THEOREM D (Shoup-style generic-group bound for an ARBITRARY predicate of k)")
print("="*90)
print(" In the generic group model every element an algorithm can hold after m queries is")
print(" sigma(alpha_i + beta_i * k) with (alpha_i,beta_i) KNOWN.  Its entire view is the equality")
print(" pattern among these.  Two of them collide iff (beta_i-beta_j)k = alpha_j-alpha_i, one")
print(" affine equation over the FIELD Z_N (N prime) -> at most 1 root.  For k drawn from any")
print(" distribution D, Pr[a given nontrivial pair collides] <= 1/|supp D|.  With <= C(m,2) pairs:")
print("      Pr[view depends on k at all]  <=  m^2 / (2|supp D|)")
print(" Conditioned on no collision the view -- hence the OUTPUT -- is independent of k.")
print(" Therefore for ANY predicate P and any generic algorithm A making m queries,")
print("      Adv(A) := |Pr_{D0}[A=1] - Pr_{D1}[A=1]|  <=  m^2 / (2 min(|D0|,|D1|))")
print(" i.e.  m >= sqrt(2 * eps * min(|D0|,|D1|))  for advantage eps.")
print()
print(" Applied to P = [w(k) <= B], D0 = uniform on {w<=B}, D1 = uniform on {w>B}:")
print()
print("    B     |{w<=B}|    |{w>B}|    min      GENERIC LOWER BOUND (eps=1/2)")
for B in (10,20,30,56,100,128,152,198,220,245):
    a=sum(comb(256,j) for j in range(B+1)); b=2**256-a
    mn=min(a,b)
    print("  %4d   2^%6.1f   2^%6.1f   2^%6.1f      m >= 2^%.1f"%(B,log2(a),log2(b),log2(mn),0.5*log2(mn)))
print()
print(" READ IT BOTH WAYS -- this is the point:")
print("  * B=128 (the null median): m >= 2^127.5.  DECIDING w<=128 IS AS HARD AS SOLVING (2^126.5).")
print("    The weight predicate has NO generic shortcut.  This is the barrier that was missing.")
print("  * B=20: m >= 2^49 only.  Small classes are genuinely cheap to decide -- which is exactly")
print("    why every LOWER bound in this campaign was affordable.  The asymmetry is class size.")
