#!/usr/bin/env python3
"""Mod-p univariate polynomial root finding (Cantor-Zassenhaus equal-degree splitting for linear
factors). Polynomials are lists of coeffs, index=degree (little-endian)."""
import random
p = 2**256 - 2**32 - 977

def norm(a):
    a=[x%p for x in a]
    while len(a)>1 and a[-1]==0: a.pop()
    return a
def pdeg(a):
    a=norm(a); return len(a)-1 if a!=[0] else -1
def padd(a,b):
    n=max(len(a),len(b)); r=[0]*n
    for i in range(len(a)): r[i]=(r[i]+a[i])%p
    for i in range(len(b)): r[i]=(r[i]+b[i])%p
    return norm(r)
def psub(a,b):
    n=max(len(a),len(b)); r=[0]*n
    for i in range(len(a)): r[i]=(r[i]+a[i])%p
    for i in range(len(b)): r[i]=(r[i]-b[i])%p
    return norm(r)
def pmul(a,b):
    a=norm(a); b=norm(b)
    if a==[0] or b==[0]: return [0]
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b):
                r[i+j]=(r[i+j]+ai*bj)%p
    return norm(r)
def pscale(a,s):
    return norm([(x*s)%p for x in a])
def pdivmod(a,b):
    a=norm(a); b=norm(b)
    db=pdeg(b)
    if db<0: raise ZeroDivisionError
    inv_lead=pow(b[-1],p-2,p)
    q=[0]*max(1,pdeg(a)-db+1); r=a[:]
    while pdeg(r)>=db and r!=[0]:
        d=pdeg(r)-db; coef=(r[-1]*inv_lead)%p
        q[d]=coef
        sub=[0]*(d)+[ (coef*bc)%p for bc in b]
        r=psub(r,sub)
    return norm(q),norm(r)
def pmod(a,b): return pdivmod(a,b)[1]
def pgcd(a,b):
    a=norm(a); b=norm(b)
    while b!=[0]:
        a,b=b,pmod(a,b)
    if a==[0]: return [0]
    return pscale(a,pow(a[-1],p-2,p))  # monic
def ppowmod(base, e, mod):
    """base^e mod (mod polynomial), e integer."""
    result=[1]; b=pmod(base,mod)
    while e>0:
        if e&1: result=pmod(pmul(result,b),mod)
        b=pmod(pmul(b,b),mod); e>>=1
    return result

def roots_mod_p(poly, tries=200):
    """Return list of distinct roots in F_p of poly (list little-endian)."""
    f=norm(poly)
    if pdeg(f)<0: return []   # zero poly: skip
    # remove zero root
    roots=[]
    # divide out leading to monic
    f=pscale(f, pow(f[-1],p-2,p))
    # g = gcd(f, x^p - x) -> product of distinct linear factors
    xp=ppowmod([0,1], p, f)      # x^p mod f
    xpx=psub(xp,[0,1])           # x^p - x
    g=pgcd(f,xpx)
    if pdeg(g)<=0:
        return []
    # recursively split g (all factors linear)
    stack=[g]; out=[]
    while stack:
        h=stack.pop()
        d=pdeg(h)
        if d<=0: continue
        if d==1:
            # h = x - r (monic) -> r = -h[0]
            out.append((-h[0])%p); continue
        # equal-degree splitting (degree 1 factors)
        split=None
        for _ in range(tries):
            r=random.randrange(p)
            b=[r,1]  # x + r
            t=ppowmod(b,(p-1)//2,h)
            t=psub(t,[1])
            c=pgcd(h,t)
            dc=pdeg(c)
            if 0<dc<d:
                split=c; break
        if split is None:
            # fallback: try to extract via gcd with x^((p-1)/2)-1 without shift
            t=psub(ppowmod([0,1],(p-1)//2,h),[1]); c=pgcd(h,t)
            if 0<pdeg(c)<d: split=c
            else:
                raise RuntimeError("split failed")
        q,_=pdivmod(h,split)
        stack.append(split); stack.append(q)
    # dedup and verify
    res=[]
    for r in set(out):
        # verify
        val=0
        for i,co in enumerate(norm(poly)):
            val=(val+co*pow(r,i,p))%p
        if val==0: res.append(r)
    return res

if __name__=='__main__':
    # self-test: (x-3)(x-5)(x-7) = x^3 -15x^2 +71x -105
    poly=pmul(pmul([(-3)%p,1],[(-5)%p,1]),[(-7)%p,1])
    print("test roots:", sorted(roots_mod_p(poly)), "expect [3,5,7]")
    # quadratic non-residue check: x^2 - a
    a=pow(random.randrange(p),2,p)
    print("sqrt test:", sorted(roots_mod_p([(-a)%p,0,1])))
