#!/usr/bin/env python3
"""Decode signed-digit HIT lines into an actual signed representation of k, and verify exactly."""
import json,sys
import numpy as np
d=json.load(open('xdata.json')); p=int(d['p']); A_=int(d['a']); N=int(d['N'])
lad=[(int(a),int(b)) for a,b in d['ladder']]
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
def nxt(s): return ((s>>1)+1)<<1
def tbl_tuples():
    """enumerate signed table tuples in build order: a=1,2,3 ; leading index even"""
    for s0 in range(0,512,2): yield (s0,)
    for s0 in range(0,512,2):
        for s1 in range(nxt(s0),512): yield (s0,s1)
    for s0 in range(0,512,2):
        for s1 in range(nxt(s0),512):
            for s2 in range(nxt(s1),512): yield (s0,s1,s2)
def val(ss):  # signed value of a tuple of signed-ladder indices
    return sum((-1 if (s&1) else 1)<<(s>>1) for s in ss)
def main():
    tf,line,tgt = sys.argv[1],sys.argv[2],(sys.argv[3] if len(sys.argv)>3 else None)
    T=(int(d['T'][0]),int(d['T'][1]))
    if tgt:
        ln=open(tgt).readline().split(); T=(int(ln[0]),int(ln[1]))
    pr=line.split(); sz=int(pr[1]); code=int(pr[2]); key=int(pr[4])
    scan=[(code>>(16*i))&0xFFFF for i in range(sz)]
    k_scan=val(scan)
    keys=np.memmap(tf,dtype=np.uint64,mode='r')
    pos=np.nonzero(keys==np.uint64(key))[0]
    print('scan indices',scan,'-> exponents/signs',[( s>>1, -1 if s&1 else 1) for s in scan])
    print('table positions matching key:',len(pos))
    tl=list(tbl_tuples())
    out=[]
    for q in pos[:50]:
        tt=tl[int(q)]
        for sgn in (1,-1):
            k = (sgn*val(tt) - k_scan) % N
            if mul(k,lad[0])==T:
                terms=[(s>>1, sgn*(-1 if s&1 else 1)) for s in tt] + [(s>>1, -(-1 if s&1 else 1)) for s in scan]
                out.append((k,sorted(terms)))
    print('VERIFIED signed representations found:',len(out))
    for k,t in out[:5]:
        print('   m=%d  terms=%s'%(len(t),t))
        print('     k mod N = %d   unsigned weight %d'%(k,bin(k).count('1')))
main()
