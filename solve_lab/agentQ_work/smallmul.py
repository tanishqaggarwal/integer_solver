#!/usr/bin/env python3
"""Q-8e: does a SMALL multiple of T land on the ladder?  m*T == 2^i*G  =>  k == 2^i * m^{-1} mod N.
Covers every k that is a ladder point divided by a small integer -- a family none of the other
sweeps reaches."""
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fastg import *
M=10**7
G=toj((31917591553801470078828036568057743875467637605644620066197178005619323650152,
       83364444556352143115103874010002344754157095926378075484791050960431190202517))
T=toj((mpz(30121525689829097248416773597728729849687459852468451992398421980273013515302),
       mpz(44859544763832475231923253825569092119321525945631045653619508440821028887)))
t0=time.time(); P=[G]
for i in range(255): P.append(jdbl(P[-1]))
lad={int(x):i for i,x in enumerate(batch_affine_x(P))}
cur=(mpz(1),mpz(1),mpz(0)); m=0; CH=1<<16
while m<M:
    buf=[]
    for _ in range(min(CH,M-m)):
        cur=jadd(cur,T); buf.append(cur)
    for j,x in enumerate(batch_affine_x(buf)):
        if x is not None and int(x) in lad:
            mm=m+j+1; i=lad[int(x)]
            k=pow(mm,-1,int(N))*pow(2,i,int(N))%int(N)
            xr=batch_affine_x([jmul(k,G)])[0]
            if xr==T[0]:
                print('*** SOLVED: m=%d i=%d  k=%d'%(mm,i,k)); sys.exit(0)
            print('x-collision m=%d i=%d (wrong sign)'%(mm,i))
    m+=len(buf)
    if m%(1<<20)==0: print('  %d  %.0fs'%(m,time.time()-t0),flush=True)
print('no small multiple of T (m <= %d) lands on the ladder.  %.0fs'%(M,time.time()-t0))
