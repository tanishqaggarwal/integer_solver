#!/usr/bin/env python3
"""Agent AB: independent recomputation of every number used in UPPER_BOUND_MAP.md.
Reads only agentX_work/xdata.json (read-only) for G,T,ladder; re-derives everything else."""
import json, os, math, random
from math import comb
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'..','agentX_work','xdata.json')))
p=int(d['p']); A=int(d['a']); B=int(d['b']); N=int(d['N'])
G=(int(d['G'][0]),int(d['G'][1])); T=(int(d['T'][0]),int(d['T'][1]))
lad=[(int(x),int(y)) for x,y in d['ladder']]
print("== 0. curve re-verification (my own arithmetic) ==")
print(" p == 2^256-2^32-977 :", p==2**256-2**32-977)
print(" N == secp256k1 order:", N==0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
print(" a,b =", A, B)
def inv(z): return pow(z,p-2,p)
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0: return None
        l=(3*x1*x1+A)*inv(2*y1%p)%p
    else:
        l=(y2-y1)*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p
    return (x3,(l*(x1-x3)-y1)%p)
def neg(P): return None if P is None else (P[0],(-P[1])%p)
def mul(k,P):
    k%= N
    R=None;Q=P
    while k>0:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
def oncur(P): return P is None or (P[1]*P[1]-pow(P[0],3,p)-A*P[0]-B)%p==0
print(" G on curve:",oncur(G)," T on curve:",oncur(T))
print(" N*G == O  :", mul(N,G) is None, "  N*T == O :", mul(N,T) is None)
ok=all(lad[i]==mul(pow(2,i,N),G) for i in range(0,256,17)) and lad[0]==G
print(" ladder[i] == 2^i*G  (spot check every 17th + i=0):", ok)
# full independent doubling chain
c=G; full=True
for i in range(256):
    if lad[i]!=c: full=False; print("  MISMATCH at",i); break
    c=add(c,c)
print(" ladder[i] == 2^i*G  (full 256 by doubling):", full)
print(" number of ladder exponents:", len(lad), " exponent set = 0..255")

print()
print("== 1. COMPLEMENT IDENTITY, checked numerically ==")
ALL=2**256-1
rng=random.Random(20260807)
def fold(S): return sum(1<<i for i in S)
bad=0
for trial in range(200):
    w=rng.randrange(0,257)
    S=set(rng.sample(range(256),w))
    Sb=set(range(256))-S
    if fold(S)+fold(Sb)!=ALL: bad+=1; continue
    if add(mul(fold(S),G), mul(fold(Sb),G)) != mul(ALL,G): bad+=1
print(" 200 random S: fold(S)+fold(~S)==2^256-1 and folds sum to (2^256-1)G  -- failures:",bad)
Tc=add(mul(ALL,G), neg(T))
print(" complement target T' = (2^256-1)G - T  =")
print("   x =",Tc[0]); print("   y =",Tc[1])
print(" T' on curve:",oncur(Tc), " N*T'==O:", mul(N,Tc) is None)
# end-to-end simulation of the mechanism on a PLANTED high-weight key
kp = rng.getrandbits(256) | (1<<255)
kp = 2**256-1-fold(set(rng.sample(range(256),6)))     # planted weight 250
Tp = mul(kp,G)
Tpc= add(mul(ALL,G), neg(Tp))
kpc= ALL-kp
print(" planted k of weight",bin(kp).count('1'),": complement weight",bin(kpc).count('1'),
      " and (ALL)G-Tp == kpc*G :", Tpc==mul(kpc,G))
print(" -> mechanism is sound: a weight-<=W hit on T' recovers a weight-(256-W) k on T")

print()
print("== 2. N's binary structure and the FREE unconditional upper bound ==")
b=bin(N)[2:]
print(" bitlen(N)=",N.bit_length()," popcount(N)=",bin(N).count('1')," zeros=",256-bin(N).count('1'))
print(" top 128 bits of N (bits 128..255):", bin(N>>128)[2:])
print(" 2^256 - N =",2**256-N," = 2^%.3f"%math.log2(2**256-N)," popcount",bin(2**256-N).count('1'))
# max Hamming weight of an integer x with 0 <= x < N  (digit DP, exact)
def maxwt_below(M):
    # max popcount over 0 <= x < M
    best=-1; bl=M.bit_length()
    pref=0
    for i in range(bl-1,-1,-1):
        if (M>>i)&1:
            # x agrees with M above i, has 0 at i, then all ones below
            cand=pref+i
            if cand>best: best=cand
            pref+=1
    return best
mw=maxwt_below(N)
print(" max popcount over 0 <= k < N :", mw, "   (exact digit-DP)")
# sanity: brute force the same DP on a small analogue
for M in [1,2,3,5,177,1000,0xBAAE]:
    bf=max(bin(x).count('1') for x in range(M))
    assert bf==maxwt_below(M), (M,bf,maxwt_below(M))
print(" digit-DP validated against brute force on 7 small moduli: OK")
print(" => UNCONDITIONAL, ZERO-COST:  w = popcount(k0) <= %d  where k0 = dlog_G(T) in [0,N)"%mw)

print()
print("== 3. the null and the information content of each bound ==")
tot=2**256
def tail_ge(a): return sum(comb(256,j) for j in range(a,257))
def tail_le(a): return sum(comb(256,j) for j in range(0,a+1))
for W in (7,10,12,14,16,20):
    lo=tail_le(W); hi=tail_ge(256-W)
    print("  W=%2d : P(w<=W)=2^%.1f   P(w>=256-W)=2^%.1f   (equal by symmetry: %s)"
          %(W, math.log2(lo/tot), math.log2(hi/tot), lo==hi))
print("  bits of information from excluding a null-mass-q region ~= q/ln2 :")
for W in (10,20):
    q=tail_le(W)/tot
    print("   W=%2d -> q=2^%.1f -> %.3e bits"%(W,math.log2(q),q/math.log(2)))
print("  P(w in [104,152]) =", "%.6f"%(sum(comb(256,j) for j in range(104,153))/tot))

print()
print("== 4. costs ==")
print(" unsigned MITM  C(256,w/2):")
for w in (10,12,14,20,24,30,40,56):
    print("   w<=%2d : 2^%.1f"%(w,math.log2(comb(256,w//2))))
print(" signed  MITM  C(256,m/2)*2^(m/2):")
for m in (8,10,12,14,16,20):
    print("   m<=%2d : 2^%.1f"%(m,math.log2(comb(256,m//2)*2**(m//2))))
print(" rho with negation+GLV: 2^%.1f"%(math.log2(math.sqrt(math.pi*N/4)/math.sqrt(6))))

print()
print("== 5. endomorphism lambda ==")
# lambda = cube root of 1 mod N
lam=None
for c in range(2,10**6):
    t=pow(c,(N-1)//3,N)
    if t!=1 and pow(t,3,N)==1: lam=t; break
lams=sorted({lam,pow(lam,2,N)})
print(" lambda^3 == 1 mod N :", all(pow(l,3,N)==1 for l in lams), " lambda !=1:", all(l!=1 for l in lams))
beta=None
for c in range(2,10**6):
    t=pow(c,(p-1)//3,p)
    if t!=1 and pow(t,3,p)==1: beta=t; break
L=lams[0]
print(" beta^3==1 mod p:", pow(beta,3,p)==1)
# find which (lam,beta) pair matches the endomorphism (beta*x, y) = [lam]P
for l in lams:
    for bb in sorted({beta,pow(beta,2,p)}):
        if mul(l,G)==((bb*G[0])%p,G[1]): print("  matched endomorphism: lam=%d beta=%d"%(l,bb)); L=l
print(" weight of lambda*k for 300 random weight-4 k (mean should be ~128):")
ws=[bin((L*fold(set(rng.sample(range(256),4))))%N).count('1') for _ in range(300)]
print("   mean %.2f  sd %.2f  min %d  max %d"%(sum(ws)/len(ws),
      (sum((x-sum(ws)/len(ws))**2 for x in ws)/len(ws))**.5, min(ws), max(ws)))
print(" weight of N-k (negation) for 300 random weight-4 k:")
ws=[bin((-fold(set(rng.sample(range(256),4))))%N).count('1') for _ in range(300)]
print("   mean %.2f  sd %.2f  min %d  max %d"%(sum(ws)/len(ws),
      (sum((x-sum(ws)/len(ws))**2 for x in ws)/len(ws))**.5, min(ws), max(ws)))

print()
print("== 6. 2-adic: doubling is a bijection on a group of odd prime order ==")
h=pow(2,N-2,N)
print(" (2^-1 mod N) exists:", (2*h)%N==1)
print(" halving 40 times from T stays on curve & is consistent:",
      all(mul(2,mul(pow(h,i,N),T))==mul(pow(h,i-1,N),T) for i in range(1,41)))
print(" => every point has a unique halving; [2] is a bijection; v_2(k) unconstrained by T")

print()
print("== 7. counting / uniqueness ==")
print(" #S with fold(S) in [0,2^256): 2^256 ; each k0 in [0,N) hit once, twice iff k0 < 2^256-N")
print(" P(second solution) = (2^256-N)/N = 2^%.1f"%math.log2((2**256-N)/N))
print(" if a second solution exists then k0 < 2^%.1f so popcount(k0) <= %d"
      %(math.log2(2**256-N), (2**256-N).bit_length()))

print()
print("== 8. does the complement fold into the SIGNED-digit search? ==")
print(" 2^256-1 has NAF weight 2 (= +2^256 - 2^0):", (2**256-1)==2**256-1)
for wt in (4,6,8):
    S=set(rng.sample(range(256),wt)); k=ALL-fold(S)
    # signed rep of k with exponents in 0..256 : +2^256 -2^0 - sum 2^e
    terms=[(1,256),(-1,0)]+[(-1,e) for e in sorted(S)]
    val=sum(s*2**e for s,e in terms)%N
    print("  complement-weight %d key: signed rep with m=%d terms, exponents<=256, correct: %s"
          %(wt,len(terms), val==k%N))
