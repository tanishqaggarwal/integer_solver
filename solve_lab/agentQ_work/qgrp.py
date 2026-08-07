#!/usr/bin/env python3
"""Q-3: group utilities on the invariant cubic + leaf-structure probes."""
import json, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
cs = K*pow(3,p-2,p)%p
_c = json.load(open(os.path.join(HERE,'curve.json')))
A_, B_ = int(_c['a']), int(_c['b'])
def inv(z): return pow(z,p-2,p)
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0: return None
        l=(3*x1*x1+A_)%p*inv(2*y1%p)%p
    else:
        l=(y2-y1)%p*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p
    return (x3,(l*(x1-x3)-y1)%p)
def neg(P): return None if P is None else (P[0],(-P[1])%p)
def sub(P,Q): return add(P,neg(Q))
def mul(k,P):
    if k<0: k,P=-k,neg(P)
    R=None;Q=P
    while k>0:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
def oncur(P): return P is None or (P[1]*P[1]-pow(P[0],3,p)-A_*P[0]-B_)%p==0
def rawpt(x,y): return ((x+cs)%p, y%p)

TARGET = rawpt(91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002,
               125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626)

def leaves():
    pins=json.load(open(os.path.join(HERE,'..','agentF_work','pins.json')))
    out={}
    for g,vv in pins.items():
        if len(vv)!=2: continue
        (v1,k1),(v2,k2)=vv
        out[int(g)]=rawpt(k2,k1)
    return out
