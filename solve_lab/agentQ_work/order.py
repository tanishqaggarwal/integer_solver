#!/usr/bin/env python3
"""Q-4: group order of the invariant cubic, computed here (Cornacchia, a=0 case), then verified
by exact scalar multiplication on random curve points."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qgrp import *
import math

def isqrt(n): return math.isqrt(n)
def sqrt_mod(a,p):
    a%=p
    if a==0: return 0
    if pow(a,(p-1)//2,p)!=1: return None
    # Tonelli-Shanks
    q=p-1; s=0
    while q%2==0: q//=2; s+=1
    if s==1: return pow(a,(p+1)//4,p)
    z=2
    while pow(z,(p-1)//2,p)!=p-1: z+=1
    m,c,t,r=s,pow(z,q,p),pow(a,q,p),pow(a,(q+1)//2,p)
    while t!=1:
        i=0; t2=t
        while t2!=1: t2=t2*t2%p; i+=1
        b=pow(c,1<<(m-i-1),p)
        m=i; c=b*b%p; t=t*c%p; r=r*b%p
    return r

# 4p = L^2 + 27 M^2
r = sqrt_mod(p-27, p)
assert r is not None, 'p !~ 1 mod 3'
x0 = r if r%2==1 else (p-r)
a_,b_ = 2*p, x0
lim = 2*isqrt(p)
while b_ > lim: a_,b_ = b_, a_%b_
L = b_; rem = 4*p - L*L
assert rem % 27 == 0
M = isqrt(rem//27); assert 27*M*M == rem, 'no rep'
print('4p = L^2 + 27 M^2 with')
print('  L =', L)
print('  M =', M)

cands = set()
for s1 in (1,-1):
    cands.add(s1*L)
    if (s1*L+9*M)%2==0: cands.add((s1*L+9*M)//2)
    if (s1*L-9*M)%2==0: cands.add((s1*L-9*M)//2)
print('candidate traces:', len(cands))

def randpt():
    while True:
        x=random.randrange(p)
        y2=(pow(x,3,p)+A_*x+B_)%p
        y=sqrt_mod(y2,p)
        if y is not None: return (x,y)
pts=[randpt() for _ in range(4)] + [TARGET]
good=[]
for t in sorted(cands):
    N=p+1-t
    if all(mul(N,P) is None for P in pts): good.append(N)
print('orders killing all test points:', good)
for N in good:
    print('  N =', N)
    print('  N == p?', N==p, ' N-p =', N-p)
    # factor small part
    n=N; f=[]
    for q in range(2,200000):
        while n%q==0: f.append(q); n//=q
        if q*q>n: break
    print('  small factors:', f, ' cofactor prime?', pow(2,n-1,n)==1 if n>1 else True, ' cofactor bits', n.bit_length())
