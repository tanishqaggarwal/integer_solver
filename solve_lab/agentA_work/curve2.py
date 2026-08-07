"""(1) Which j=0 curve is y^2 = x^3 + b1 (b1 from the PINNED point (x12186,x16742))?
   (2) Do the instance's large literals pair into points on it?
   (3) Where do my two congruence residues sit?"""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
N_SEC=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
x1=v[12186]%P; y1=v[16742]%P
b1=(y1*y1-pow(x1,3,P))%P
print('x1=%d\ny1=%d\nb1=%d'%(x1,y1,b1))
print('point on y^2=x^3+b1 :', (y1*y1-pow(x1,3,P)-b1)%P==0)
# is b1/7 a sixth power?  (=> curve isomorphic to secp256k1)
def is_kth_power(a,k):
    a%=P
    if a==0: return True
    g=(P-1)//__import__('math').gcd(P-1,k)
    return pow(a,(P-1)//__import__('math').gcd(P-1,k),P)==1
import math
for k in [2,3,6]:
    e=(P-1)//math.gcd(P-1,k)
    print('b1/7 is a %d-th power mod p: %s'%(k, pow(b1*pow(7,P-2,P)%P, e, P)==1))
# curve arithmetic on y^2 = x^3 + b1
def add(Pt,Qt):
    if Pt is None: return Qt
    if Qt is None: return Pt
    x1_,y1_=Pt; x2_,y2_=Qt
    if x1_==x2_ and (y1_+y2_)%P==0: return None
    if Pt==Qt: lam=3*x1_*x1_%P*pow(2*y1_,P-2,P)%P
    else: lam=(y2_-y1_)*pow(x2_-x1_,P-2,P)%P
    x3=(lam*lam-x1_-x2_)%P; y3=(lam*(x1_-x3)-y1_)%P
    return (x3,y3)
def mul(k,Pt):
    R=None; Q=Pt
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
G=(x1,y1)
print('N_secp * (x1,y1) =', mul(N_SEC,G))
# the six possible orders for j=0 curves: p+1-t, t from CM
# find a,b with 4p = a^2+27b^2 by brute search over b using Cornacchia-lite
def cornacchia(d,n):
    # solve x^2 + d y^2 = n
    for r in [pow(-d,(n+1)//4,n)] if n%4==3 else []:
        pass
    return None
orders=set()
# derive orders empirically: order must satisfy m*G = O for the true order m; test divisors
cands=[P+1, N_SEC]
print('checking candidate orders by m*G:')
for m in cands:
    print('   m=%d -> %s'%(m, mul(m,G)))
# 6 twists orders sum relations: use the fact that #E_b depends only on b's class mod 6th powers
# brute: find the order by testing (p+1-t) for t from the CM equation 4p=a^2+27c^2
import itertools
found=None
a2=None
# solve 4p = a^2 + 27 c^2
fourp=4*P
c=math.isqrt(fourp//27)
while c>0:
    rem=fourp-27*c*c
    a=math.isqrt(rem)
    if a*a==rem:
        found=(a,c); break
    c-=1
print('4p = a^2+27c^2 with (a,c) =', found)
if found:
    a,c=found
    ts=set()
    for s1 in (1,-1):
        for s2 in (1,-1):
            ts.add(s1*a); ts.add(s1*(a+9*c)//2 if (a+9*c)%2==0 else None)
            ts.add(s1*(a-9*c)//2 if (a-9*c)%2==0 else None)
    ts={t for t in ts if t is not None}
    for t in sorted(ts):
        m=P+1-t
        print('   order candidate %d (t=%d): m*G = %s  prime? %s'%(m,t,mul(m,G),
              'yes' if __import__('sympy').isprime(m) else 'no'))
