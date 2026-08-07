#!/usr/bin/env python3
"""Q-8c: all k whose ON-bits lie inside a window of W=34 consecutive ladder positions.
k = a*2^s, a < 2^W.  MITM: a = alo + 2^17*ahi."""
import sys,os,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fastg import *
import numpy as np
W=34; H=17
TX=mpz(30121525689829097248416773597728729849687459852468451992398421980273013515302)
G=toj((31917591553801470078828036568057743875467637605644620066197178005619323650152,
       83364444556352143115103874010002344754157095926378075484791050960431190202517))
T=toj((TX,mpz(44859544763832475231923253825569092119321525945631045653619508440821028887)))
t0=time.time(); P=[G]
for i in range(255): P.append(jdbl(P[-1]))
SS=256-W+1
n=SS*(1<<H)
def build(neg,shift0,tag):
    keys=np.empty(n,dtype=np.uint64); codes=np.empty(n,dtype=np.int64); m=0
    for s in range(SS):
        base=P[s+shift0] if s+shift0<256 else None
        buf=[];cb=[]
        if base is None: continue
        acc=(mpz(1),mpz(1),mpz(0))
        Bs=jneg(base) if neg else base
        cur=T if neg else (mpz(1),mpz(1),mpz(0))
        for a in range(1<<H):
            buf.append(cur); cb.append(s*(1<<H)+a)
            cur=jadd(cur,Bs)
            if len(buf)==1<<16:
                for j,x in enumerate(batch_affine_x(buf)): keys[m+j]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
                codes[m:m+len(buf)]=cb; m+=len(buf); buf=[];cb=[]
        if buf:
            for j,x in enumerate(batch_affine_x(buf)): keys[m+j]=(int(x)&0xFFFFFFFFFFFFFFFF) if x is not None else 0
            codes[m:m+len(buf)]=cb; m+=len(buf)
        if s%40==0: print('  %s s=%d  %.0fs'%(tag,s,time.time()-t0),flush=True)
    return keys[:m],codes[:m]
kA,cA=build(False,0,'A'); oA=np.argsort(kA); sA=kA[oA]
print('A done %.0fs'%(time.time()-t0),flush=True)
kB,cB=build(True,H,'B')
print('B done %.0fs'%(time.time()-t0),flush=True)
pos=np.searchsorted(sA,kB); tried=0
for m in range(len(kB)):
    pp=pos[m]
    while pp<len(sA) and sA[pp]==kB[m]:
        ca=int(cA[oA[pp]]); cb=int(cB[m])
        sa,alo=divmod(ca,1<<H); sb,ahi=divmod(cb,1<<H)
        if sa==sb:
            k=(alo+(ahi<<H))<<sa
            if 0<k<(1<<256):
                x=batch_affine_x([jmul(k,G)])[0]
                if x==TX: print('*** SOLVED: k =',k); sys.exit(0)
        tried+=1; pp+=1
print('no k confined to a %d-bit window.  collisions %d  %.0fs'%(W,tried,time.time()-t0))
