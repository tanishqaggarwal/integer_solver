from fractions import Fraction as Fr
import json
p=2**256-2**32-977

# ---- exact-integer LLL (size reduction + Lovasz swaps, Gram-Schmidt over Fractions) ----
def lll(B, delta=Fr(3,4)):
    B=[row[:] for row in B]; n=len(B)
    def dot(u,v): return sum(a*b for a,b in zip(u,v))
    def gram():
        Bs=[]; mu=[[Fr(0)]*n for _ in range(n)]
        for i in range(n):
            bi=[Fr(x) for x in B[i]]
            for j in range(i):
                mu[i][j]=Fr(dot(B[i],[int(x) for x in [0]*0] ) ) if False else dot([Fr(x) for x in B[i]],Bs[j])/dot(Bs[j],Bs[j])
                bi=[bi[k]-mu[i][j]*Bs[j][k] for k in range(len(bi))]
            Bs.append(bi)
        return Bs,mu
    Bs,mu=gram(); k=1
    while k<n:
        for j in range(k-1,-1,-1):
            if abs(mu[k][j])>Fr(1,2):
                q=int((mu[k][j]+Fr(1,2)).__floor__()) if mu[k][j]>=0 else -int((-mu[k][j]+Fr(1,2)).__floor__())
                B[k]=[B[k][t]-q*B[j][t] for t in range(len(B[k]))]
                Bs,mu=gram()
        if dot(Bs[k],Bs[k])>=(delta-mu[k][k-1]**2)*dot(Bs[k-1],Bs[k-1]):
            k+=1
        else:
            B[k],B[k-1]=B[k-1],B[k]; Bs,mu=gram(); k=max(k-1,1)
    return B

# Genuine lattice problem: the gadget's exact target x_31731 must be 0, on the G1=G2=0 coset.
# x_9118 = x7068 - (7376877*p)*s ,  x_8731 = x4432 - p*t   (s,t integers = x_17325, x_9413)
# x_31731(s,t) = c31 + s31*s + t31*t  (from _crtsolve). Want minimize |x_31731|.
# This is a 1-eq integer minimization -> lattice {(s,t): } find (s,t) with x_31731 near 0.
G=json.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/gadget_affine.json'))
b={int(k):int(v) for k,v in G['b'].items()}
a9118={int(k):int(v) for k,v in G['s9118'].items()}
a8731={int(k):int(v) for k,v in G['s8731'].items()}
x7068=int(G['x7068']); x4432=int(G['x4432'])
M1=7376877*p; M2=p
def form(q):
    const=b[q]+a9118[q]*x7068+a8731[q]*x4432
    return const, -a9118[q]*M1, -a8731[q]*M2
c31,s31,t31=form(31731)
# Build lattice to minimize |c31 + s31*s + t31*t| with s,t moderate: embed
# rows: [s31*W, s, 0], [t31*W, 0, t-scale], target row [c31*W,0,0]
W=1  # weight on the target coordinate
# Use CVP-style: basis {(s31, 1, 0),(t31,0,1)}; target (-c31,0,0). Reduce [ [s31,1,0],[t31,0,1],[c31,0,0] ] won't embed target as lattice pt.
# Standard: reduce B=[[s31,K,0],[t31,0,K]] then we want lattice pt near (-c31,0,0): do embedding
K=1
B=[[s31, K, 0],[t31, 0, K],[c31, 0, 0]]  # last row = the shift; LLL will combine
R=lll(B)
print("LLL-reduced basis (gadget exact-target lattice):")
best=None
# search small combinations of reduced rows to minimize |first coord| (=x_31731) with the c31 offset
import itertools
rows=R
for coeffs in itertools.product(range(-3,4),repeat=3):
    val=[sum(coeffs[i]*rows[i][j] for i in range(3)) for j in range(3)]
    # we need the combination to represent c31 + s31*s+t31*t form: the '1' in col1/col2 track s,t, col0 tracks x_31731
    # Only combos where col using the c31 row exactly once represent a real (s,t)
    pass
# Simpler: directly report gcd wall
from math import gcd
g=gcd(s31,t31)
print(f"gcd(s31,t31)={g}")
print(f"x_31731 achievable values = c31 + (multiples of g). c31 mod g = {c31%g}")
print(f"min |x_31731| on G1=G2=0 coset = {min(c31%g, g-(c31%g))}  (0 achievable? {c31%g==0})")
print(f"=> x_31731=0 EXACTLY is {'ACHIEVABLE' if c31%g==0 else 'IMPOSSIBLE'} while G1=G2=0")
# Also the mod-p residue wall for x_2239:
c22,s22,t22=form(2239)
print(f"x_2239 mod p on coset: coefficients s,t ≡ {s22%p},{t22%p} mod p (both 0 => residue FIXED)")
print(f"x_2239 residue = {c22%p}  (needs 0 for gadget) => {'OK' if c22%p==0 else 'WALL'}")
