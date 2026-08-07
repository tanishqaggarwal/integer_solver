"""Extract A,B and the six labelled coordinates from a frame's residual, and map the
three points onto secp256k1."""
import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
from sympy.ntheory.residue_ntheory import nthroot_mod
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
inv3=pow(3,-1,P)
Bs = 64019533680030876408443198762210829058751700634554282185987325820393598524794
n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
Gp=(Gx,Gy)
_u=nthroot_mod((Bs*pow(7,-1,P))%P,6,P)
IU2=pow(_u*_u%P,-1,P); IU3=pow(pow(_u,3,P),-1,P)
def tosec(x,y): return (((x+K*inv3)%P)*IU2%P, y*IU3%P)
def add(A,B):
    if A is None: return B
    if B is None: return A
    (xa,ya),(xb,yb)=A,B
    if xa==xb and (ya+yb)%P==0: return None
    l=(3*xa*xa*pow(2*ya,-1,P)%P) if A==B else ((yb-ya)*pow((xb-xa)%P,-1,P)%P)
    xc=(l*l-xa-xb)%P; return (xc,(l*(xa-xc)-ya)%P)
def neg(A): return None if A is None else (A[0],(-A[1])%P)
def sub(A,B): return add(A,neg(B))
def mul(k,A):
    R=None; k%=n
    while k:
        if k&1: R=add(R,A)
        A=add(A,A); k>>=1
    return R
c11,c12,c21,c22=8646263,1073965,10159099,6926539
_det=(c11*c22-c12*c21)%P; _inv=pow(_det,-1,P)

def pencil(f1,f2):
    """given symbolic a19297,a19299 return (A,B) polynomials"""
    def lc(a,b):
        out={}
        for m,c in f1.items(): out[m]=(out.get(m,0)+c*a)%P
        for m,c in f2.items(): out[m]=(out.get(m,0)+c*b)%P
        return {m:c for m,c in out.items() if c}
    return lc(c22*_inv%P,(-c12)*_inv%P), lc((-c21)*_inv%P,c11*_inv%P)

def label(Apoly,Bpoly,NB):
    """identify x1,y1,x2,y2,x3,y3 from the monomial structure"""
    varsA={NB[k] for m in Apoly for k,_ in m}
    varsB={NB[k] for m in Bpoly for k,_ in m}
    cubed={NB[k] for m in Apoly for k,e in m if e==3}
    deg3mons=[m for m in Apoly if sum(e for _,e in m)==3]
    in3={NB[k] for m in deg3mons for k,_ in m}
    x3set=in3-cubed
    y3set=varsB-varsA
    yset=varsA-in3
    if len(cubed)!=2 or len(x3set)!=1 or len(y3set)!=1 or len(yset)!=2: return None
    x3=x3set.pop(); y3=y3set.pop()
    ix={u:i for i,u in enumerate(NB)}
    # sign of the y3*X monomial: +1 -> x2, -1 -> x1  (B = y3(x2-x1)+...)
    x1=x2=None
    for m,c in Bpoly.items():
        ks=[NB[k] for k,e in m for _ in range(e)]
        if y3 in ks and len(ks)==2:
            other=ks[0] if ks[1]==y3 else ks[1]
            if c%P==1: x2=other
            elif (-c)%P==1: x1=other
    if x1 is None or x2 is None: return None
    # y1 pairs with x2 (coefficient +1), y2 pairs with x1 (coefficient -1)
    y1=y2=None
    for m,c in Bpoly.items():
        ks=[NB[k] for k,e in m for _ in range(e)]
        if len(ks)!=2: continue
        s=set(ks)
        if x2 in s and (s-{x2}) <= yset and c%P==1: y1=(s-{x2}).pop()
        if x1 in s and (s-{x1}) <= yset and (-c)%P==1: y2=(s-{x1}).pop()
    if y1 is None or y2 is None: return None
    return dict(x1=x1,y1=y1,x2=x2,y2=y2,x3=x3,y3=y3)
