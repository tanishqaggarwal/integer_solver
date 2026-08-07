"""Fast EC over y^2=x^3+B (a=0), Jacobian coords + Montgomery batch inversion."""
import numpy as np
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
B=64019533680030876408443198762210829058751700634554282185987325820393598524794
N=115792089237316195423570985008687907852837564279074904382605163141518161494337  # prime group order

def jdouble(P):
    X,Y,Z=P
    if Y==0: return (0,0,0)
    A=X*X%p; Bq=Y*Y%p; C=Bq*Bq%p
    D=2*(((X+Bq)*(X+Bq)-A-C)%p)%p
    E=3*A%p; F=E*E%p
    X3=(F-2*D)%p; Y3=(E*(D-X3)-8*C)%p; Z3=2*Y*Z%p
    return (X3,Y3,Z3)

def jadd_affine(P,q):
    """P Jacobian, q=(x2,y2) affine."""
    X1,Y1,Z1=P
    if Z1==0: return (q[0],q[1],1)
    ZZ=Z1*Z1%p
    U2=q[0]*ZZ%p
    S2=q[1]*Z1%p*ZZ%p
    H=(U2-X1)%p
    r=(S2-Y1)%p
    if H==0:
        if r==0: return jdouble(P)
        return (0,0,0)
    HH=H*H%p; I=4*HH%p; J=H*I%p; r2=2*r%p; V=X1*I%p
    X3=(r2*r2-J-2*V)%p
    Y3=(r2*(V-X3)-2*Y1*J)%p
    Z3=((Z1+H)*(Z1+H)-ZZ-HH)%p
    return (X3,Y3,Z3)

def batch_norm(js):
    """list of Jacobian -> list of affine (x,y) ; skips Z==0 as None."""
    zs=[t[2] for t in js]
    pref=[1]*(len(zs)+1)
    acc=1
    for i,z in enumerate(zs):
        pref[i]=acc
        if z: acc=acc*z%p
    inv=pow(acc,p-2,p)
    out=[None]*len(js)
    for i in range(len(zs)-1,-1,-1):
        z=zs[i]
        if not z: continue
        zi=inv*pref[i]%p
        inv=inv*z%p
        zi2=zi*zi%p
        out[i]=(js[i][0]*zi2%p, js[i][1]*zi2%p*zi%p)
    return out

# affine helpers
def aadd(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if x1==x2:
        if (y1+y2)%p==0: return None
        l=3*x1*x1%p*pow(2*y1%p,p-2,p)%p
    else:
        l=(y2-y1)*pow((x2-x1)%p,p-2,p)%p
    x3=(l*l-x1-x2)%p
    return (x3,(l*(x1-x3)-y1)%p)
def aneg(P): return None if P is None else (P[0],(-P[1])%p)
def amul(k,P):
    k%=N
    R=None; Q=P
    while k:
        if k&1: R=aadd(R,Q)
        Q=aadd(Q,Q); k>>=1
    return R
MASK=(1<<62)-1
def key(x): return np.int64(x & MASK)
