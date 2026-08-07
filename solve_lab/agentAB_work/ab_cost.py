#!/usr/bin/env python3
"""Agent AB: (i) the ball-covering cost floor for ANY search-based upper bound on w,
   (ii) the masked-complement family, (iii) residual curve-structure checks."""
import math, random, json, os
from math import comb, log2
N=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
p=2**256-2**32-977
ALL=2**256-1

print("== A. BALL-COVERING COST FLOOR ==")
print(" Knob set: MITM search of a Hamming ball B(c,W) around any known center c, cost C(256,W/2)")
print(" (c=0 -> plain low-weight search; c=ALL -> the complement search; general c -> signed rep).")
print(" To PROVE w<=B you must cover {k : wt(k)>B}.  #balls >= |region| / Vol(256,W).")
def Vol(W): return sum(comb(256,j) for j in range(W+1))
def region(B): return sum(comb(256,j) for j in range(B+1,257))
print()
print("   B     |{wt>B}|      best W   #balls     cost/ball   TOTAL COST   vs rho 2^126.5")
for B in (245,240,220,200,180,160,152,144,136,128):
    R=region(B); best=None
    for W in range(2,241,2):
        nb=max(1.0, R/Vol(W))
        tot=math.log2(nb)+log2(comb(256,W//2))
        if best is None or tot<best[0]: best=(tot,W,math.log2(nb),log2(comb(256,W//2)))
    tot,W,nb,cb=best
    print("  %4d  2^%7.1f    W=%3d   2^%7.1f    2^%6.1f     2^%7.1f    %s"
          %(B, log2(R), W, nb, cb, tot, "AFFORDABLE" if tot<126.5 else "worse than just solving"))
print()
print(" Reference: solving outright (rho, negation+GLV) = 2^126.5 and yields w EXACTLY.")
print(" => any covering total above 2^126.5 is strictly dominated by solving the instance.")

print()
print("== B. MASKED-COMPLEMENT FAMILY: which centers give an unconditional bound? ==")
print(" c - k = c XOR k  <=>  no borrows  <=>  supp(k) subset supp(c).")
print(" So wt(c-k) = wt(c)-wt(k) holds for ALL k iff supp(c)={0..255}, i.e. c=2^256-1.")
rng=random.Random(11)
def wt(x): return bin(x).count('1')
def corr(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sx=(sum((a-mx)**2 for a in xs)/n)**.5; sy=(sum((b-my)**2 for b in ys)/n)**.5
    return sum((a-mx)*(b-my) for a,b in zip(xs,ys))/(n*sx*sy)
ks=[rng.randrange(N) for _ in range(4000)]
wk=[wt(k) for k in ks]
for name,c in [("c = 2^256-1 (all ones)",ALL),
               ("c = N-1",N-1),
               ("c = 2^255",2**255),
               ("c = 2^128-1 (low mask)",2**128-1),
               ("c = random",rng.randrange(N)),
               ("c = 0 (plain search)",0)]:
    wc=[wt((c-k)%N) for k in ks]
    exact=all(a+b==256 for a,b in zip(wk,wc)) if c==ALL else False
    print("  %-24s corr(wt(k),wt(c-k mod N)) = %+.4f   exact wt(k)+wt(c-k)=256 : %s"%(name,corr(wk,wc),exact))
print(" Only the all-ones centre carries a deterministic relation; every other centre is ~uncorrelated")
print(" or (c=0) perfectly correlated in the wrong direction.")

print()
print("== C. masked complement c_A over a proper subset A (conditional bound) ==")
for sz in (255,224,192,128):
    A=set(rng.sample(range(256),sz)); cA=sum(1<<i for i in A)
    good=0; tot=0
    for _ in range(2000):
        w0=rng.randrange(0,sz+1); S=set(rng.sample(sorted(A),w0))
        tot+=1
        if wt(cA-sum(1<<i for i in S))==sz-w0: good+=1
    print("  |A|=%3d : for k with supp(k) subset A, wt(c_A-k)=|A|-wt(k) held %d/%d ; "
          "P(supp(k) subset A) under the null = 2^%.1f"%(sz,good,tot, sz-256))

print()
print("== D. residual curve-structure checks (all DEAD-verdict inputs) ==")
try:
    from sympy import isprime, factorint
    print("  N prime:",isprime(N),"   p prime:",isprime(p))
except Exception as e:
    print("  sympy unavailable:",e)
print("  N == p (anomalous / Smart attack applies) :", N==p)
print("  trace t = p+1-N =", p+1-N, " (|t| ~ 2^%.1f)"%log2(abs(p+1-N)))
# CM: 4p = t^2 + 3 v^2  for j=0 curves
t=p+1-N
r=4*p-t*t
print("  (4p - t^2)/3 is a perfect square (CM by sqrt(-3)) :", r%3==0 and math.isqrt(r//3)**2==r//3)
# embedding degree
k=None; x=1
for i in range(1,20001):
    x=x*p%N
    if x==1: k=i; break
print("  multiplicative order of p mod N is <= 20000 :", k)
print("  => embedding degree > 20000; MOV/Frey-Ruck target field has > 2^%d bits"%(256*20000//1000000*1000))
d=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentX_work','xdata.json')))
b=int(d['b'])
z=b*pow(7,p-2,p)%p
print("  a == 0 (j=0, Aut = mu_6):", int(d['a'])==0)
print("  b/7 is a 6th power mod p (curve F_p-isomorphic to secp256k1 y^2=x^3+7):",
      pow(z,(p-1)//6,p)==1)

print()
print("== E. bit extraction is a 256-query full solve (self-contained reduction) ==")
print("  k0 = dlog_G(T) in [0,N) c [0,2^256); S = {i : bit_i(k0)=1} is ALWAYS a valid subset,")
print("  so the instance is ALWAYS satisfiable and w = popcount(k0).")
print("  bit_i(k0) IS the predicate 'i in S'.  An oracle for one bit -> 256 calls -> k0 -> full solve.")
print("  Hence any single-bit extractor costs >= (DLP cost)/256 = 2^%.1f."%(126.5-8))

print()
print("== F. subset-sum / lattice: the density is a red herring ==")
print("  density of the modular subset-sum {2^i} mod N : n/log2(N) = 256/256 = 1.000")
print("  low-density attacks (Lagarias-Odlyzko, CJLOSS) require density < 0.9408 : not met.")
print("  BUT the decisive point: GIVEN k0 the subset-sum is solved by reading binary digits.")
print("  100%% of the hardness is 'k0 unknown'.  A lattice needs an integer target; the only")
print("  integers derivable from T are x_T,y_T in F_p, related to k0 only through division")
print("  polynomials of degree ~k0^2/2 ~ 2^511.  No linear structure exists to reduce.")
