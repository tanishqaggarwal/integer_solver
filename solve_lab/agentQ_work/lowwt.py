#!/usr/bin/env python3
"""Q-8b: meet-in-the-middle for a LOW-HAMMING-WEIGHT k.  P_i = 2^i G, i=0..255.
Side A = sums of <=3 ladder points; side B = T - (sums of <=3).  Collision => weight<=6 k."""
import sys, os, time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fastg import *
import numpy as np
TX=mpz(30121525689829097248416773597728729849687459852468451992398421980273013515302)
G=toj((31917591553801470078828036568057743875467637605644620066197178005619323650152,
       83364444556352143115103874010002344754157095926378075484791050960431190202517))
T=toj((TX,mpz(44859544763832475231923253825569092119321525945631045653619508440821028887)))
t0=time.time()
P=[G]
for i in range(255): P.append(jdbl(P[-1]))
NPTS=1+256+256*255//2+256*255*254//6
print('ladder built; per side %d points  %.0fs'%(NPTS,time.time()-t0),flush=True)

def gen(base,sign):
    Q=[jneg(q) for q in P] if sign<0 else P
    yield base,-1
    for i in range(256): yield jadd(base,Q[i]), i
    for i in range(256):
        Bi=jadd(base,Q[i])
        for j in range(i+1,256): yield jadd(Bi,Q[j]), 1000+i*256+j
    for i in range(256):
        Bi=jadd(base,Q[i])
        for j in range(i+1,256):
            Bij=jadd(Bi,Q[j])
            for k in range(j+1,256): yield jadd(Bij,Q[k]), 1000000+(i*65536+j*256+k)

def build(base,sign,tag):
    keys=np.empty(NPTS,dtype=np.uint64); codes=np.empty(NPTS,dtype=np.int64)
    buf=[]; cb=[]; n=0
    for pt,c in gen(base,sign):
        buf.append(pt); cb.append(c)
        if len(buf)==1<<16:
            for j,x in enumerate(batch_affine_x(buf)): keys[n+j]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
            codes[n:n+len(buf)]=cb; n+=len(buf); buf=[]; cb=[]
            if n%(1<<20)<(1<<16): print('  %s %d  %.0fs'%(tag,n,time.time()-t0),flush=True)
    if buf:
        for j,x in enumerate(batch_affine_x(buf)): keys[n+j]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
        codes[n:n+len(buf)]=cb; n+=len(buf)
    return keys[:n],codes[:n]

def decode(c):
    c=int(c)
    if c<0: return []
    if c<1000: return [c]
    if c<1000000: c-=1000; return [c//256,c%256]
    c-=1000000; return [c//65536,(c//256)%256,c%256]

kA,cA=build((mpz(1),mpz(1),mpz(0)),1,'A')
oA=np.argsort(kA); sA=kA[oA]
print('side A sorted %.0fs'%(time.time()-t0),flush=True)
kB,cB=build(T,-1,'B')
print('side B built %.0fs'%(time.time()-t0),flush=True)
pos=np.searchsorted(sA,kB); hits=0
for m in range(len(kB)):
    pp=pos[m]
    while pp<len(sA) and sA[pp]==kB[m]:
        a=decode(cA[oA[pp]]); b=decode(cB[m]); S=set(a)|set(b)
        if len(S)==len(a)+len(b) and S:
            k=sum(1<<i for i in S)
            x=batch_affine_x([jmul(k,G)])[0]
            if x is not None and x==TX:
                print('*** SOLVED: k =',k,' bits',sorted(S)); sys.exit(0)
        hits+=1; pp+=1
print('no k of Hamming weight <= 6.  collisions inspected:',hits,'  %.0fs'%(time.time()-t0))
