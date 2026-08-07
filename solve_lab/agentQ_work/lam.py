#!/usr/bin/env python3
"""Q-8f: the cubic has a=0, so it carries an efficient endomorphism phi(X,Y)=(beta*X,Y) with
beta^3=1 mod p, acting as multiplication by lambda (lambda^3=1 mod N).  Two tests:
 (1) is T = +-lambda^j * 2^i * G  (1536 candidates)?
 (2) MITM for k = a + b*lambda mod N with |a|,|b| < 2^21."""
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fastg import *
import numpy as np
TX=mpz(30121525689829097248416773597728729849687459852468451992398421980273013515302)
G=toj((31917591553801470078828036568057743875467637605644620066197178005619323650152,
       83364444556352143115103874010002344754157095926378075484791050960431190202517))
T=toj((TX,mpz(44859544763832475231923253825569092119321525945631045653619508440821028887)))
# cube roots of unity
def cuberoots(mod):
    r=[]
    g=2
    while True:
        c=pow(g,(int(mod)-1)//3,int(mod))
        if c!=1: break
        g+=1
    return [1,c,c*c%int(mod)]
BET=cuberoots(p); LAM=cuberoots(N)
beta=mpz(BET[1]); lam=LAM[1]
gx=batch_affine_x([G])[0]
phiG=(beta*G[0]%p, G[1], G[2])
ok = batch_affine_x([jmul(lam,G)])[0]==batch_affine_x([phiG])[0]
if not ok: lam=LAM[2]; ok = batch_affine_x([jmul(lam,G)])[0]==batch_affine_x([phiG])[0]
print('endomorphism phi(X,Y)=(beta X, Y) acts as multiplication by lambda:',ok)
t0=time.time(); P=[G]
for i in range(255): P.append(jdbl(P[-1]))
# (1)
for i in range(256):
    for j in range(3):
        for s in (1,-1):
            k=(s*pow(2,i,int(N))*pow(lam,j,int(N)))%int(N)
            x=batch_affine_x([jmul(k,G)])[0]
            if x==TX: print('*** SOLVED (1): k =',k); sys.exit(0)
print('(1) T is not +-lambda^j * 2^i * G')
# (2)
Wb=1<<21
LG=jmul(lam,G)
keys=np.empty(2*Wb,dtype=np.uint64); codes=np.empty(2*Wb,dtype=np.int64)
m=0; cur=(mpz(1),mpz(1),mpz(0)); buf=[];cb=[]
for a in range(Wb):
    buf.append(cur); cb.append(a)
    buf.append(jneg(cur)); cb.append(-a)
    cur=jadd(cur,G)
    if len(buf)>=(1<<16):
        for j2,x in enumerate(batch_affine_x(buf)): keys[m+j2]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
        codes[m:m+len(buf)]=cb; m+=len(buf); buf=[];cb=[]
for j2,x in enumerate(batch_affine_x(buf)): keys[m+j2]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
codes[m:m+len(buf)]=cb; m+=len(buf)
keys=keys[:m]; codes=codes[:m]; o=np.argsort(keys); sk=keys[o]
print('side aG built %.0fs'%(time.time()-t0),flush=True)
cur=T; buf=[];cb=[]
def flush(buf,cb):
    if not buf: return
    xs=batch_affine_x(buf)
    arr=np.array([(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0 for x in xs],dtype=np.uint64)
    pos=np.searchsorted(sk,arr)
    for t in range(len(arr)):
        pp=pos[t]
        while pp<len(sk) and sk[pp]==arr[t]:
            k=(int(codes[o[pp]])+cb[t]*lam)%int(N)
            if k and batch_affine_x([jmul(k,G)])[0]==TX:
                print('*** SOLVED (2): k =',k); sys.exit(0)
            pp+=1
for b in range(Wb):
    buf.append(cur); cb.append(b)
    buf.append(jadd(T,jneg(jadd(cur,jneg(T))))); cb.append(-b)
    cur=jadd(cur,jneg(LG))
    if len(buf)>=(1<<16): flush(buf,cb); buf=[];cb=[]
flush(buf,cb)
print('(2) no k = a + b*lambda with |a|,|b| < 2^21.  %.0fs'%(time.time()-t0))
