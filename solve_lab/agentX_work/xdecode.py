#!/usr/bin/env python3
"""Decode a HIT line into the actual subset S, and verify it exactly."""
import json,os,sys,itertools
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'xdata.json')))
p=int(d['p']); A_=int(d['a']); N=int(d['N'])
lad=[(int(a),int(b)) for a,b in d['ladder']]
Treal=(int(d['T'][0]),int(d['T'][1]))
def inv(z):return pow(z,p-2,p)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x1,y1=P;x2,y2=Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0:return None
        l=(3*x1*x1+A_)%p*inv(2*y1%p)%p
    else: l=(y2-y1)%p*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
def mul(k,P):
    R=None;Q=P
    while k>0:
        if k&1:R=add(R,Q)
        Q=add(Q,Q);k>>=1
    return R
def unrank(pos):
    """table file position -> subset (sizes 1..4, lexicographic blocks)"""
    B=[256,32640,2763520,174792640]
    for s,(n) in enumerate(B,start=1):
        if pos<n:
            # unrank the pos-th s-combination of range(256) in lex order
            res=[];start=0;rem=pos
            from math import comb
            for t in range(s):
                for v in range(start,256):
                    c=comb(255-v,s-t-1)
                    if rem<c: res.append(v);start=v+1;break
                    rem-=c
            return res
        pos-=n
    raise ValueError
def decode_scan(code,size):
    return [(code>>(8*i))&0xFF for i in range(size)]
def main():
    tblfile=sys.argv[1]; line=sys.argv[2]; Tx=sys.argv[3] if len(sys.argv)>3 else None
    parts=line.split()
    size=int(parts[1]); code=int(parts[2]); key=int(parts[4])
    beta=decode_scan(code,size)
    k=np.memmap(tblfile,dtype=np.uint64,mode='r')
    hits=np.nonzero(k==np.uint64(key))[0]
    print('scan subset beta =',beta)
    print('table positions with that key:',hits[:10],'count',len(hits))
    T=Treal
    if Tx:
        dd=json.load(open(Tx)); T=(int(dd['T'][0]),int(dd['T'][1]))
    for h in hits[:20]:
        alpha=unrank(int(h))
        kk=sum(1<<i for i in alpha)+sum(1<<i for i in beta)
        P=mul(kk%N,lad[0])
        ok = P==T
        print('  alpha=%s  k=%d  bits(k)=%s  weight=%d  k<2^256:%s  k*G==T:%s'%(
            alpha,kk,sorted(i for i in range(300) if (kk>>i)&1),bin(kk).count('1'),kk<2**256,ok))
main()
