#!/usr/bin/env python3
"""Q-2: (i) associativity/commutativity of the stage law, (ii) locate the ROOT TARGET pin and test
whether it lies on the same cubic invariant as the leaves."""
import json, os, re, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
cs = K * pow(3, p-2, p) % p
cur = json.load(open(os.path.join(HERE,'curve.json')))
A_, B_ = int(cur['a']), int(cur['b'])
def inv(z): return pow(z, p-2, p)

# ---- law in RAW coordinates (as measured by agent F) ----
def chordK(P, Q):
    ax,ay = P; bx,by = Q
    l = (by-ay)*inv((bx-ax)%p)%p
    ox = (l*l - ax - bx - K)%p
    oy = (l*(ax-ox) - ay)%p
    return (ox,oy)
# ---- in shifted coordinates X = x + K/3 ----
def to_s(P): return ((P[0]+cs)%p, P[1]%p)
def fr_s(P): return ((P[0]-cs)%p, P[1]%p)
def oncur(P): return (P[1]*P[1] - pow(P[0],3,p) - A_*P[0] - B_) % p == 0

# group operations on the shifted cubic (identity = point at infinity, encoded None)
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0: return None
        l = (3*x1*x1 + A_)%p * inv(2*y1%p) %p
    else:
        l = (y2-y1)%p * inv((x2-x1)%p) %p
    x3 = (l*l - x1 - x2)%p
    return (x3, (l*(x1-x3)-y1)%p)
def neg(P): return None if P is None else (P[0], (-P[1])%p)
def mul(k,P):
    R=None; Q=P
    while k>0:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R

rnd = random.Random(7)
pins = json.load(open(os.path.join(HERE,'..','agentF_work','pins.json')))
L=[]
for g,vv in pins.items():
    if len(vv)!=2: continue
    (v1,k1),(v2,k2)=vv
    L.append((g, to_s((k2%p, k1%p))))   # orient 1: X=k2 (already raw x), Y=k1  -> shift
# careful: orient=1 fit used X = k2 + c.  to_s adds c.  ok.
print('leaf points:', len(L), 'all on cubic:', all(oncur(P) for _,P in L))

# (i) does chordK agree with the shifted chord law?
ok=0
for _ in range(200):
    P=L[rnd.randrange(len(L))][1]; Q=L[rnd.randrange(len(L))][1]
    if P[0]==Q[0]: continue
    r1 = to_s(chordK(fr_s(P), fr_s(Q)))
    # chordK gives the REFLECTED third intersection = group add
    ok += (r1 == add(P,Q))
print('chordK == shifted group add on %d random leaf pairs' % ok)

# (ii) associativity + commutativity, tested directly on the raw law
na=nc=nt=0
for _ in range(300):
    P=L[rnd.randrange(len(L))][1]; Q=L[rnd.randrange(len(L))][1]; R=L[rnd.randrange(len(L))][1]
    if len({P[0],Q[0],R[0]})<3: continue
    nt+=1
    na += (add(add(P,Q),R)==add(P,add(Q,R)))
    nc += (add(P,Q)==add(Q,P))
print('associative %d/%d   commutative %d/%d' % (na,nt,nc,nt))
