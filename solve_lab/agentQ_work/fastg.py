#!/usr/bin/env python3
"""Fast group arithmetic on the invariant cubic (Jacobian coords + gmpy2 + batch inversion)."""
import gmpy2
from gmpy2 import mpz
p  = mpz(115792089237316195423570985008687907853269984665640564039457584007908834671663)
A_ = mpz(0)
B_ = mpz(64019533680030876408443198762210829058751700634554282185987325820393598524794)
N  = mpz(115792089237316195423570985008687907852837564279074904382605163141518161494337)
K  = mpz(97553848499418123410591666447050222001188385549510401465815187079080512838891)
CS = K*gmpy2.invert(mpz(3),p) % p          # shift  X = x + K/3

def jdbl(P):
    X,Y,Z=P
    if Y==0: return (mpz(1),mpz(1),mpz(0))
    S=(4*X*Y*Y)%p; M=(3*X*X)%p
    X3=(M*M-2*S)%p
    return (X3,(M*(S-X3)-8*pow(Y,4,p))%p,(2*Y*Z)%p)
def jadd(P,Q):
    X1,Y1,Z1=P; X2,Y2,Z2=Q
    if Z1==0: return Q
    if Z2==0: return P
    Z1Z1=Z1*Z1%p; Z2Z2=Z2*Z2%p
    U1=X1*Z2Z2%p; U2=X2*Z1Z1%p
    S1=Y1*Z2*Z2Z2%p; S2=Y2*Z1*Z1Z1%p
    H=(U2-U1)%p; r=(S2-S1)%p
    if H==0:
        if r==0: return jdbl(P)
        return (mpz(1),mpz(1),mpz(0))
    HH=H*H%p; HHH=H*HH%p; V=U1*HH%p
    X3=(r*r-HHH-2*V)%p
    return (X3,(r*(V-X3)-S1*HHH)%p,(Z1*Z2*H)%p)
def jmul(k,P):
    R=(mpz(1),mpz(1),mpz(0))
    k=int(k)
    while k>0:
        if k&1: R=jadd(R,P)
        P=jdbl(P); k>>=1
    return R
def jneg(P): return (P[0],(-P[1])%p,P[2])
def toj(P): return (mpz(P[0]),mpz(P[1]),mpz(1)) if P else (mpz(1),mpz(1),mpz(0))

def batch_affine_x(pts):
    """affine x-coordinates of a list of Jacobian points (Z=0 -> None)."""
    zs=[P[2] for P in pts]
    pre=[]; acc=mpz(1)
    for z in zs:
        pre.append(acc)
        if z: acc=acc*z%p
    inv=gmpy2.invert(acc,p)
    out=[None]*len(pts)
    for i in range(len(pts)-1,-1,-1):
        z=zs[i]
        if not z: continue
        zi=inv*pre[i]%p
        inv=inv*z%p
        out[i]=pts[i][0]*zi*zi%p
    return out
