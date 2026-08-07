#!/usr/bin/env python3
"""Q-8a: baby-step giant-step for a SMALL discrete log  k*G = T  (and for -T, i.e. k near N).
Covers k < 2^BABY+GIANT.  A hit fully solves the instance."""
import sys, os, time, json
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fastg import *
import numpy as np
BABY=1<<22; GIANT=1<<22
G_raw=(31917591553801470078828036568057743875467637605644620066197178005619323650152,
       83364444556352143115103874010002344754157095926378075484791050960431190202517)
T_raw=(30121525689829097248416773597728729849687459852468451992398421980273013515302,
       44859544763832475231923253825569092119321525945631045653619508440821028887)
G=toj(G_raw); T=toj(T_raw)
t0=time.time()
print('building %d baby steps...'%BABY, flush=True)
keys=np.empty(BABY,dtype=np.uint64)
CH=1<<16
P=(mpz(1),mpz(1),mpz(0))
idx=0
while idx<BABY:
    n=min(CH,BABY-idx); buf=[]
    for _ in range(n):
        P=jadd(P,G); buf.append(P)
    xs=batch_affine_x(buf)
    for j,x in enumerate(xs):
        keys[idx+j]= (int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
    idx+=n
    if idx%(1<<20)==0: print('  %d  %.0fs'%(idx,time.time()-t0),flush=True)
order=np.argsort(keys); skeys=keys[order]
print('baby table done %.0fs'%(time.time()-t0),flush=True)
S=jmul(BABY,G)
for sign,TT in (('+',T),('-',jneg(T))):
    Q=TT; found=None
    i=0
    while i<GIANT:
        n=min(CH,GIANT-i); buf=[]
        for _ in range(n):
            buf.append(Q); Q=jadd(Q,jneg(S))
        xs=batch_affine_x(buf)
        arr=np.array([(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0 for x in xs],dtype=np.uint64)
        pos=np.searchsorted(skeys,arr)
        for j in range(n):
            pp=pos[j]
            while pp<len(skeys) and skeys[pp]==arr[j]:
                b=int(order[pp])+1
                k=(i+j)*BABY+b
                if sign=='-': k=int(N)-k
                R=jmul(k,G)
                xr=batch_affine_x([R])[0]
                if xr is not None and xr==mpz(T_raw[0]):
                    found=k; break
                pp+=1
            if found: break
        if found: break
        i+=n
        if i%(1<<20)==0: print('  giant %s %d  %.0fs'%(sign,i,time.time()-t0),flush=True)
    if found:
        print('*** DLP SOLVED: k =',found)
        json.dump({'k':str(found)},open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'k_found.json'),'w'))
        sys.exit(0)
print('no small k in the searched range; %.0fs'%(time.time()-t0))
