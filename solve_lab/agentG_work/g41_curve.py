"""Re-test the curve hypothesis in GENERAL Weierstrass form, deriving every constant
from my own A,B identities (nothing imported).

My exact identities (g22_numbers.py, derived from the symbolic model):
   A = (x1-x2)^2 * (x3 + x1 + x2 + K) - (y2-y1)^2
   B = y3*(x2-x1) + y1*x2 - x1*y2 + x3*y2 - x3*y1     [ = (y3+y1)(x2-x1) - (x1-x3)(y2-y1) ]
On  y^2 = x^3 + a2 x^2 + a4 x + a6  the chord-and-tangent addition is
   lam = (y2-y1)/(x2-x1),  x3 = lam^2 - a2 - x1 - x2,  y3 = lam(x1-x3) - y1
i.e. exactly  A = 0 with K = a2, and B = 0.  So K IS the a2 coefficient.
a4, a6 then follow from requiring P1 and P2 to lie on the curve (2 linear equations).
j = 0  <=>  a4 = a2^2/3.
"""
import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v=L.load(src); ad.fwd(v,rounds=6)
C={'x1':22649,'y1':16742,'x2':14853,'y2':31339,'x3':22162,'y3':30213}
val={k:v[u]%P for k,u in C.items()}
for k in ['x1','y1','x2','y2','x3','y3']: print('%s = x%-6d = %d'%(k,C[k],val[k]))
x1,y1,x2,y2,x3,y3=[val[k] for k in ['x1','y1','x2','y2','x3','y3']]
print('\nK = a2 =',K)
# a4, a6 from P1, P2 on  y^2 = x^3 + K x^2 + a4 x + a6
r1=(y1*y1 - x1**3 - K*x1*x1)%P
r2=(y2*y2 - x2**3 - K*x2*x2)%P
# r1 = a4*x1 + a6 ; r2 = a4*x2 + a6
a4=((r1-r2)*pow((x1-x2)%P,-1,P))%P
a6=(r1-a4*x1)%P
print('a4 =',a4)
print('a6 =',a6)
print('\n--- j = 0 test:  a4 == a2^2/3 ? ---')
inv3=pow(3,-1,P)
t=(K*K*inv3)%P
print('a2^2/3 =',t)
print('MATCH (j=0):', a4==t)
print('\n--- short form after x -> x - a2/3 ---')
Ashort=(a4 - K*K*inv3)%P
Bshort=(a6 - K*a4*inv3 + 2*pow(K,3,P)*pow(27,-1,P))%P
print('A_short =',Ashort)
print('B_short =',Bshort,'   (secp256k1 has A=0, B=7)')
print('B_short == 7 ?', Bshort==7)
# is B_short a 6th power (=> isomorphic to secp256k1 over F_p)?
if Ashort==0 and Bshort:
    q=(P-1)//6
    u=(Bshort*pow(7,-1,P))%P
    print('(B_short/7) is a 6th power ?', pow(u,q,P)==1)
    print('(B_short/7) is a square    ?', pow(u,(P-1)//2,P)==1)
    print('(B_short/7) is a cube      ?', pow(u,(P-1)//3,P)==1)
print('\n--- is P3 on the same curve? ---')
print('P3 residual =',(y3*y3 - x3**3 - K*x3*x3 - a4*x3 - a6)%P)
print('\n--- group order test on the SHORT curve y^2 = x^3 + A_short x + B_short ---')
n_secp=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
def add(Pt,Q,a):
    if Pt is None: return Q
    if Q is None: return Pt
    (xa,ya),(xb,yb)=Pt,Q
    if xa==xb and (ya+yb)%P==0: return None
    if Pt==Q: l=(3*xa*xa+a)*pow(2*ya,-1,P)%P
    else: l=(yb-ya)*pow((xb-xa)%P,-1,P)%P
    xc=(l*l-xa-xb)%P; yc=(l*(xa-xc)-ya)%P
    return (xc,yc)
def mul(k,Pt,a):
    R=None
    while k:
        if k&1: R=add(R,Pt,a)
        Pt=add(Pt,Pt,a); k>>=1
    return R
sh=lambda x:(x+K*inv3)%P   # x_short = x + a2/3
for nm,(xx,yy) in [('P1',(x1,y1)),('P2',(x2,y2))]:
    Q=(sh(xx),yy)
    res=(yy*yy - pow(Q[0],3,P) - Ashort*Q[0] - Bshort)%P
    print('%s on short curve? residual=%d'%(nm,res))
    if res==0:
        print('   [n_secp]%s = %s'%(nm, 'O (INFINITY)' if mul(n_secp,Q,Ashort) is None else mul(n_secp,Q,Ashort)))
