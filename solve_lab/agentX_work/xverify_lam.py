#!/usr/bin/env python3
"""Agent X: independent re-verification of the endomorphism facts and the orbit exclusion."""
import json
d=json.load(open('xdata.json')); p=int(d['p']); A_=int(d['a']); B_=int(d['b']); N=int(d['N'])
lad=[(int(a),int(b)) for a,b in d['ladder']]; T=(int(d['T'][0]),int(d['T'][1]))
def inv(z):return pow(z,p-2,p)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x1,y1=P;x2,y2=Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0:return None
        l=(3*x1*x1+A_)%p*inv(2*y1%p)%p
    else: l=(y2-y1)%p*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p;return (x3,(l*(x1-x3)-y1)%p)
def mul(k,P):
    k%=N; R=None;Q=P
    while k>0:
        if k&1:R=add(R,Q)
        Q=add(Q,Q);k>>=1
    return R
# cube roots of unity
betas=[b for b in (pow(g,(p-1)//3,p) for g in (2,3,5,7,11)) if b!=1]
beta=betas[0]; assert pow(beta,3,p)==1 and beta!=1
lams=[l for l in (pow(g,(N-1)//3,N) for g in (2,3,5,7,11)) if l!=1]
G=lad[0]
lam=None
for l in set(lams):
    if mul(l,G)==(beta*G[0]%p,G[1]): lam=l;break
    if mul(l,G)==(pow(beta,2,p)*G[0]%p,G[1]): lam=l;beta=pow(beta,2,p);break
print('endomorphism phi(X,Y)=(beta*X,Y) == multiplication by lambda :', lam is not None)
print('  lambda^3 == 1 mod N :', pow(lam,3,N)==1)
# check on 5 random-ish ladder points
ok=all(mul(lam,lad[i])==(beta*lad[i][0]%p, lad[i][1]) for i in (0,17,64,199,255))
print('  verified on 5 ladder points :', ok)
# orbit exclusion: is T = +-lambda^j * 2^i * G ?
orb=set()
for i in range(256):
    x,y=lad[i]
    for j in range(3):
        xx=pow(beta,j,p)*x%p
        orb.add((xx,y)); orb.add((xx,(-y)%p))
print('T in the 1536-point orbit {+- lambda^j 2^i G} :', T in orb, ' (orbit size %d)'%len(orb))
print('sqrt(3) speedup only: |<lambda>| = 3, so the endomorphism divides the search space by 3, i.e. 2^128 -> 2^127.2')
