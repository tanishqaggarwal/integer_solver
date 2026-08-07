#!/usr/bin/env python3
"""Q-8d: MITM for Hamming weight <= 7:  side A = sums of <=3 ladder pts, side B = T - sums of <=4."""
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fastg import *
import numpy as np
TX=mpz(30121525689829097248416773597728729849687459852468451992398421980273013515302)
G=toj((31917591553801470078828036568057743875467637605644620066197178005619323650152,
       83364444556352143115103874010002344754157095926378075484791050960431190202517))
T=toj((TX,mpz(44859544763832475231923253825569092119321525945631045653619508440821028887)))
t0=time.time(); P=[G]
for i in range(255): P.append(jdbl(P[-1]))
def genA():
    yield (mpz(1),mpz(1),mpz(0)),0
    for i in range(256): yield P[i],1<<i
    for i in range(256):
        Bi=P[i]
        for j in range(i+1,256):
            Bij=jadd(Bi,P[j]); yield Bij,(1<<i)|(1<<j)
            for k in range(j+1,256): yield jadd(Bij,P[k]),(1<<i)|(1<<j)|(1<<k)
NA=1+256+256*255//2+256*255*254//6
keys=np.empty(NA,dtype=np.uint64); codes=[]
m=0; buf=[]
for pt,c in genA():
    buf.append(pt); codes.append(c)
    if len(buf)==1<<16:
        for j,x in enumerate(batch_affine_x(buf)): keys[m+j]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
        m+=len(buf); buf=[]
        if m%(1<<21)<(1<<16): print('  A %d %.0fs'%(m,time.time()-t0),flush=True)
for j,x in enumerate(batch_affine_x(buf)): keys[m+j]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
m+=len(buf); keys=keys[:m]
oA=np.argsort(keys); sA=keys[oA]
print('A built %d  %.0fs'%(m,time.time()-t0),flush=True)
Q=[jneg(q) for q in P]
cnt=0
buf=[];cb=[]
def flush(buf,cb):
    global cnt
    if not buf: return
    xs=batch_affine_x(buf)
    arr=np.array([(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0 for x in xs],dtype=np.uint64)
    pos=np.searchsorted(sA,arr)
    for t in range(len(arr)):
        pp=pos[t]
        while pp<len(sA) and sA[pp]==arr[t]:
            ka=int(codes[oA[pp]]); kk=ka|cb[t]
            if ka & cb[t]==0 and kk:
                x=batch_affine_x([jmul(kk,G)])[0]
                if x==TX: print('*** SOLVED: k =',kk); sys.exit(0)
            pp+=1
    cnt+=len(buf)
for i in range(256):
    Bi=jadd(T,Q[i])
    buf.append(Bi); cb.append(1<<i)
    for j in range(i+1,256):
        Bij=jadd(Bi,Q[j]); buf.append(Bij); cb.append((1<<i)|(1<<j))
        for k in range(j+1,256):
            Bijk=jadd(Bij,Q[k]); buf.append(Bijk); cb.append((1<<i)|(1<<j)|(1<<k))
            for l in range(k+1,256):
                buf.append(jadd(Bijk,Q[l])); cb.append((1<<i)|(1<<j)|(1<<k)|(1<<l))
                if len(buf)>=(1<<16): flush(buf,cb); buf=[];cb=[]
    if i%8==0: print('  B i=%d done=%d %.0fs'%(i,cnt,time.time()-t0),flush=True)
flush(buf,cb)
print('no k of Hamming weight <= 7.  %.0fs'%(time.time()-t0))
