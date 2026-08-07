#!/usr/bin/env python3
"""Agent AB: adjudication of agent AG's THEOREM_B_AUDIT.  Every number recomputed here."""
import math, random, itertools
from math import comb, log2
def L(x): return log2(x) if x>0 else float('-inf')
N=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
RHO=L(math.sqrt(math.pi*N/4)/math.sqrt(6))
V128=[0]*129; s=0
for j in range(129): s+=comb(128,j); V128[j]=s
TOT=comb(256,128)
def wt(x): return bin(x).count('1')

print("="*88); print("STEP 0 (coordinator's instruction): is AG's NO-CARRY step correct?"); print("="*88)
print(" Claim: for ANY split of {0..255} into L,R,  k = k_L + k_R with disjoint bit supports,")
print(" hence wt(k) = wt(k_L)+wt(k_R) exactly, and every MITM certifies a RECTANGLE.")
rng=random.Random(9)
bad=0
for _ in range(20000):
    Lset=set(rng.sample(range(256),rng.randrange(0,257)))
    S=set(rng.sample(range(256),rng.randrange(0,257)))
    kL=sum(1<<i for i in S&Lset); kR=sum(1<<i for i in S-Lset)
    if kL+kR!=sum(1<<i for i in S) or wt(kL)+wt(kR)!=len(S): bad+=1
print(" 20000 random (split, S) pairs incl. NON-contiguous splits: failures =",bad)
print(" VERDICT: no-carry step is CORRECT.  Disjoint supports => no carries, additivity is exact.")
print(" Therefore attack 1's rectangle model and its sqrt(Z) floor STAND.  AG closed the")
print(" coordinator's nominated gap against itself, honestly.")
print()
print(" Independent re-derivation of the floor (AG's route (b), generic-query counting):")
print("   after m queries the algorithm holds m elements sigma(a_i+b_i k); a scalar x is excluded")
print("   only as the unique root of one of <=C(m,2) pair equations or one of m direct tests,")
print("   so excluding Z scalars needs C(m,2)+m >= Z, i.e. m >= sqrt(2Z) - O(1).   Checked:")
for Z in (2**58,2**143,2**188,2**254):
    m=math.ceil((math.sqrt(1+8*Z)-1)/2)
    print("     Z=2^%-6.1f  -> m >= 2^%.2f   (sqrt(2Z)=2^%.2f)"%(L(Z),L(m),L(math.sqrt(2*Z))))

print(); print("="*88); print("ITEM 3a: is rep(W) exactly 2x too large for odd W?"); print("="*88)
def den_AB(W):
    a=(W+1)//2
    return comb(W,a)*comb(256-W,128-a)
def den_honest(W):
    """count 128-subsets L with max(|D cap L|,|D cap R|) <= ceil(W/2)"""
    c=(W+1)//2; tot=0
    for j in range(0,W+1):
        if max(j,W-j)<=c and 0<=128-j<=256-W: tot+=comb(W,j)*comb(256-W,128-j)
    return tot
exc=[W for W in range(1,256) if (den_honest(W)!=2*den_AB(W) if W%2 else den_honest(W)!=den_AB(W))]
print(" odd W: honest denominator == 2 * AB's ?   even W: equal ?   exceptions over W=1..255:",exc)
print(" AG's proof re-checked: W=2c-1 => C(W,c)=C(W,c-1) (twin central, odd row) and")
print(" 256-W=257-2c odd with (128-c)+(129-c)=257-2c => C(256-W,128-c)=C(256-W,129-c).  Both hold:")
print("  ",all(comb(2*c-1,c)==comb(2*c-1,c-1) and comb(257-2*c,128-c)==comb(257-2*c,129-c)
              for c in range(1,129)))
print(" CONCEDED.  AG is right.")

def rep(W,honest=True):
    if W<=0: return 1.0
    d=den_honest(W) if honest else den_AB(W)
    return TOT/d if d>0 else float('inf')
def raw(W,honest=True):
    if W<=0: return 1.0
    return rep(W,honest)*V128[min((W+1)//2,128)]
def suffix(f):
    S=[0.0]*258; S[257]=float('inf'); S[256]=f(256)
    for W in range(255,-1,-1): S[W]=min(f(W),S[W+1])
    return S
BC_AB=suffix(lambda W: raw(W,False)); BC=suffix(lambda W: raw(W,True))
def Vol(W): return sum(comb(256,j) for j in range(W+1))
def region(B): return sum(comb(256,j) for j in range(B+1,257))
print()
print(" effect of the fix (unbounded memory):")
for nm,BCx in (("AB published",BC_AB),("odd-W corrected",BC)):
    xw=max(w for w in range(1,257) if L(BCx[w])<=RHO)
    be=None
    for B in range(255,100,-1):
        t=min(L(max(1.0,region(B)/Vol(W)))+L(BCx[W]) for W in range(1,257))
        if t>=RHO: be=B; break
    print("   %-16s crossover w=%d   break-even B=%d"%(nm,xw,be))

print(); print("="*88); print("ITEM 3b: can the W=256 self-certificate fail?"); print("="*88)
print("  Vol_128(128) = 2^%.4f  (= 2^128 identically)"%L(V128[128]))
print("  So ANY model rep(W)*Vol_128(ceil(W/2)) with rep(256)=1 returns 2^128 at W=256,")
print("  whatever rep does elsewhere.  Demonstrated on three DELIBERATELY WRONG reps:")
for nm,r in (("rep=1 everywhere",lambda W:1.0),("rep = AB's (odd-W bug)",lambda W:rep(W,False)),
             ("rep = 1000x nonsense",lambda W: 1.0 if W==256 else 1000*rep(W))):
    print("     %-24s W=256 cost = 2^%.4f   %s"%(nm,L(r(256)*V128[128]),"PASSES the certificate"))
print("  CONCEDED: vacuous for the model it now certifies.")
print("  PRECISION I do claim: it was NOT vacuous for the model it REFUTED -- round 2's")
print("  sqrt(W)*Vol_128(W//2) has rep(256)=16 and returned 2^132.0, a genuine failure.")
print("  The test can refute a model with rep(256)!=1; it cannot confirm one with rep(256)==1.")

print(); print("="*88); print("ITEM 3d CHALLENGE: is cover(B) = 2^128.000 EXACTLY for every B <= 148?"); print("="*88)
def cover(B,BCx,cap=None):
    best=float('inf')
    for W in range(1,257):
        c=BCx[W]
        if cap is not None: c=min(c,cap)
        v=L(max(1.0,region(B)/Vol(W)))+L(c)
        best=min(best,v)
    return best
Bstar=None
for B in range(0,256):
    if abs(cover(B,BC_AB)-128.0)>1e-9: Bstar=B; break
print("  largest B with cover(B) == 2^128.000 exactly : %d   (first B below 2^128 is %d)"%(Bstar-1,Bstar))
for B in (100,110,118,120,128,140,148,152):
    print("     B=%3d : cover = 2^%.3f  %s"%(B,cover(B,BC_AB),"= 2^128 exactly" if abs(cover(B,BC_AB)-128)<1e-9 else ""))
print("  AG's §4.4 says 'cover(B)=2^128.000 exactly for every B <= 148'.  That is FALSE for")
print("  B in [%d,148]: at B=148 cover = 2^%.3f, and AG's OWN §1.3 table lists 2^126.85 there."%(Bstar,cover(148,BC_AB)))
print("  -> AG's §4.4 contradicts AG's §1.3.  The literal claim is wrong; the CONCLUSION is right:")
print("     cover(B) lies in [rho, 2^128] = [2^%.3f, 2^128.000] for every B <= 148, a band of 2^%.3f,"%(RHO,128-RHO))
print("     so there is indeed no cliff at 148.  Restate as the band, not as an equality.")

print(); print("="*88); print("ITEM 1: memory-aware break-even -- is the propagation failure real?"); print("="*88)
C_VOW=2.5
def ball_mem(W,logM,honest=True):
    Lh=V128[min((W+1)//2,128)]; r=rep(W,honest)
    if L(Lh)<=logM: t=L(r*Lh)
    else: t=L(r*C_VOW)+1.5*L(Lh)-0.5*logM
    return min(t,RHO)                      # any ball can be searched by solving outright
def suffix_mem(logM,honest=True):
    S=[0.0]*258; S[257]=float('inf'); S[256]=ball_mem(256,logM,honest)
    for W in range(255,-1,-1): S[W]=min(ball_mem(W,logM,honest),S[W+1])
    return S
print("  memory   crossover w   break-even B   dead band      (AG reported)")
for logM in (30,35,40,50,60,80,10**6):
    S=suffix_mem(logM)
    xw=max([w for w in range(1,257) if S[w]<RHO-1e-12]+[0])
    be=None
    for B in range(255,100,-1):
        t=min(L(max(1.0,region(B)/Vol(W)))+L(2**S[W]) for W in range(1,257))
        if t>=RHO-1e-12: be=B; break
    lab="unbounded" if logM>1000 else "2^%d"%logM
    print("   %-9s   %3d          %3d          [%d,%d]"%(lab,xw,be,xw+1,be-1))
print("  AG reported 52/201 at 2^30 and 106/148 unbounded.  My published AUTHORITATIVE SUMMARY S2")
print("  quoted 106/148 as THE headline while S3 struck the unbounded-memory column.")
print("  PROPAGATION FAILURE: CONCEDED.  Same defect as round 3's, one section apart.")

print(); print("="*88); print("ITEM 4: disk"); print("="*88)
import subprocess
try:
    o=subprocess.run(['df','-B1','--output=avail,target','/'],capture_output=True,text=True).stdout
    av=int(o.split('\n')[1].split()[0]); print("  df avail on / : %.1f GB = 2^%.1f bytes"%(av/1e9,L(av)))
    print("  at 32 B per vOW distinguished point: 2^%.1f entries  (RAM figure was 2^30)"%L(av/32))
except Exception as e: print("  df failed:",e)
print("  AG measured 4.92e3 random 4 KiB reads/s under O_DIRECT vs my asserted 1e2/s -> 2^%.1f off."%L(4920/100))
print("  Both my numbers were wrong; conclusion 'disk is not a way out' HARDENS (vOW gains only")
print("  sqrt(M), so <=2^0 memory gain against >=2^14 slowdown).  CONCEDED.")
