"""Map my three points onto secp256k1 proper and identify them."""
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
src='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v=L.load(src); ad.fwd(v,rounds=6)
C={'P1':(22649,16742),'P2':(14853,31339),'P3':(22162,30213)}
t=(Bs*pow(7,-1,P))%P
u=nthroot_mod(t,6,P)
print('u =',u,' check u^6==B/7:',pow(u,6,P)==t)
iu2=pow(u*u%P,-1,P); iu3=pow(pow(u,3,P),-1,P)
def tosec(x,y): return ((x+K*inv3)%P*iu2%P, y*iu3%P)
def add(A,B,a=0):
    if A is None: return B
    if B is None: return A
    (xa,ya),(xb,yb)=A,B
    if xa==xb and (ya+yb)%P==0: return None
    if A==B: l=3*xa*xa*pow(2*ya,-1,P)%P
    else: l=(yb-ya)*pow((xb-xa)%P,-1,P)%P
    xc=(l*l-xa-xb)%P; yc=(l*(xa-xc)-ya)%P
    return (xc,yc)
def mul(k,A):
    R=None; k%=n
    while k:
        if k&1: R=add(R,A)
        A=add(A,A); k>>=1
    return R
def neg(A): return None if A is None else (A[0],(-A[1])%P)
pts={}
Gp=(Gx,Gy)
for nm,(ix,iy) in C.items():
    Q=tosec(v[ix]%P, v[iy]%P)
    pts[nm]=Q
    ok=(Q[1]*Q[1]-pow(Q[0],3,P)-7)%P==0
    print('%s -> secp256k1 (%d, %d)  oncurve=%s  [n]Q=%s'%(nm,Q[0],Q[1],ok,'O' if mul(n,Q) is None else 'NOT O'))
print('\nG =',Gp)
print('P1 == G ?',pts['P1']==Gp, ' P2 == G ?',pts['P2']==Gp, ' P3 == G ?',pts['P3']==Gp)
print('P1 == -G ?',pts['P1']==neg(Gp))
S=add(pts['P1'],pts['P2'])
print('\nP1+P2 =',S)
print('P1+P2 == P3 ?',S==pts['P3'])
print('P3-P1 == P2 ?', add(pts['P3'],neg(pts['P1']))==pts['P2'])
print('\n--- small-multiple identification (k up to 4096) ---')
R=None
tab={}
for k in range(1,4097):
    R=add(R,Gp); tab[R]=k
for nm,Q in pts.items():
    hits=[]
    if Q in tab: hits.append(('%s = [%d]G'%(nm,tab[Q])))
    if neg(Q) in tab: hits.append('%s = -[%d]G'%(nm,tab[neg(Q)]))
    print(nm, hits if hits else 'not a small multiple of G')
D=add(pts['P3'],neg(add(pts['P1'],pts['P2'])))
print('\nP3 - (P1+P2) =',D)
print('  is it a small multiple of G?', tab.get(D), tab.get(neg(D)))
