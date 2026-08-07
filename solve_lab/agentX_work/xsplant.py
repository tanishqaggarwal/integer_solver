#!/usr/bin/env python3
"""Plant a SIGNED-DIGIT target: k = sum eps_j 2^{e_j}. argv: 'e:sign,e:sign,...' out.txt"""
import json,sys
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
terms=[]
for t in sys.argv[1].split(','):
    e,s=t.split(':'); terms.append((int(e),int(s)))
k=sum(s*(1<<e) for e,s in terms)
P=mul(k%N,lad[0])
open(sys.argv[2],'w').write('%d %d\n'%P + ''.join('%d %d\n'%(x,y) for x,y in lad))
print('planted terms=%s  m=%d'%(terms,len(terms)))
print('  k mod N = %d'%(k%N))
print('  UNSIGNED Hamming weight of (k mod N) = %d   <-- what the plain sweep would need'%bin(k%N).count('1'))
# expected scan-side s0 for a given split: scan part = the top b terms
if len(sys.argv)>3:
    b=int(sys.argv[3]); scan=sorted(terms)[-b:]
    e0,eps0=scan[0]
    s0=2*e0+(1 if eps0>0 else 0)   # sigma = -eps ; s=2e if sigma=+1 (eps=-1), 2e+1 if sigma=-1 (eps=+1)
    print('  split: table=%s  scan=%s  -> expected scan s0 = %d'%(sorted(terms)[:-b],scan,s0))
    print('S0=%d'%s0)
