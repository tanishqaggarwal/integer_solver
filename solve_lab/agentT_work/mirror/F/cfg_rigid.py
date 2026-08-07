#!/usr/bin/env python3
"""Configuration-conditional mod-p rigidity.

Fix a BOOLEAN configuration (which of the boolean free inputs are 1).  Wires whose free-input support
lies inside the boolean free inputs then have values determined by that configuration alone, so the
bilinear selector products become LINEAR with known coefficients.  With that, build:
  Z_cfg : wires == 0 (mod p) for EVERY assignment with this boolean configuration
  a weighted union-find over wires from every definition / residual atom that reduces to alpha*u = beta*v
Everything derived is valid for every assignment with the given boolean configuration.
"""
import sys,os,re,json,pickle,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
SUP=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))

# boolean-constrained free inputs
BOOLFREE=set()
for a in E.res:
    m=re.fullmatch(r'\(\(x(\d+)\*x(\d+)\)-x(\d+)\)',a) or re.fullmatch(r'\(x(\d+)\*\(1-x(\d+)\)\)',a) \
      or re.fullmatch(r'\(x(\d+)\*\(x(\d+)-1\)\)',a)
    if m and len(set(m.groups()))==1: BOOLFREE.add(int(m.group(1)))
FREE=set(E.free)
BOOLFREE&=FREE
BCONE=[v for v in range(NV) if set(SUP[str(v)])<=BOOLFREE]
BCONESET=set(BCONE)

defrhs={}
for a in E.order:
    k,v,rhs=E.cls[a]; defrhs[v]=rhs

def build(on_bits):
    """on_bits: iterable of boolean free inputs set to 1 (all other free inputs arbitrary)."""
    v=[0]*NV
    for b in on_bits: v[b]=1
    E.run(v)                                   # values of BCONE wires are now correct & config-only
    val={w:v[w] for w in BCONE}
    # ---- Z closure (== 0 mod p for every assignment with this boolean config) ----
    Z=set(w for w in BCONE if val[w]%p==0)
    def isZ(n):
        o=n[0]
        if o=='c': return n[1]%p==0
        if o=='v': return n[1] in Z
        if o=='neg': return isZ(n[1])
        if o=='*': return isZ(n[1]) or isZ(n[2])
        return isZ(n[1]) and isZ(n[2])
    changed=True
    while changed:
        changed=False
        for w,rhs in defrhs.items():
            if w in Z: continue
            if isZ(rhs): Z.add(w); changed=True
    # ---- linearisation ----
    def lin(n):
        """(coef, var) if n == coef*x_var (mod p) for every assignment with this config; var None => 0."""
        o=n[0]
        if o=='v':
            w=n[1]
            if w in Z: return (0,None)
            if w in BCONESET: return None if val[w]%p else (0,None)   # a known nonzero constant, not a wire
            return (1,w)
        if o=='c': return (0,None) if n[1]%p==0 else None
        if o=='neg':
            r=lin(n[1]); return None if r is None else ((-r[0])%p,r[1])
        if o=='*':
            L,R=n[1],n[2]
            for X,Y in ((L,R),(R,L)):
                # X a known constant (BCONE wire or literal), Y linear
                c=None
                if X[0]=='c': c=X[1]%p
                elif X[0]=='v' and X[1] in BCONESET: c=val[X[1]]%p
                if c is not None:
                    if c==0: return (0,None)
                    r=lin(Y)
                    return None if r is None else ((c*r[0])%p, r[1])
            a=lin(L); b=lin(R)
            if a is not None and (a[1] is None or a[0]==0): return (0,None)
            if b is not None and (b[1] is None or b[0]==0): return (0,None)
            return None
        a=lin(n[1]); b=lin(n[2])
        if a is None or b is None: return None
        sg=1 if o=='+' else -1
        if a[1] is None or a[0]==0: return ((sg*b[0])%p, b[1])
        if b[1] is None or b[0]==0: return (a[0]%p, a[1])
        if a[1]==b[1]:
            c=(a[0]+sg*b[0])%p
            return (c,a[1]) if c else (0,None)
        return None
    # ---- weighted union-find ----
    par=list(range(NV)); fac=[1]*NV
    def find(x):
        root=x; f=1
        while par[root]!=root: f=f*fac[root]%p; root=par[root]
        # path compress
        cur=x; cf=f
        while par[cur]!=cur:
            nxt=par[cur]; nf=cf*pow(fac[cur],p-2,p)%p
            par[cur]=root; fac[cur]=cf
            cur=nxt; cf=nf
        return root,f
    def union(a,b,c):
        ra,fa=find(a); rb,fb=find(b)
        if ra==rb: return (fa*pow(fb,p-2,p))%p==c%p
        par[ra]=rb; fac[ra]=c%p*fb%p*pow(fa,p-2,p)%p
        return True
    conf=0; nl=0
    for w,rhs in defrhs.items():
        if w in Z or w in BCONESET: continue
        r=lin(rhs)
        if r is None or r[1] is None: continue
        if not union(w,r[1],r[0]): conf+=1
        nl+=1
    forced={}     # wire -> value mod p forced to a constant
    nr=0
    for a in E.res:
        ast=E.atoms[a]
        if ast[0] not in ('-','+'): continue
        L,R=ast[1],ast[2]
        rl=lin(R)
        if rl is None or not (rl[1] is None or rl[0]==0): continue    # R must vanish mod p
        # L == 0 mod p
        if L[0] in ('-','+'):
            a1=lin(L[1]); a2=lin(L[2]); sg=1 if L[0]=='+' else -1
            if a1 is None or a2 is None: continue
            if a1[1] is not None and a2[1] is not None and a1[0] and a2[0]:
                c=(-sg*a2[0]*pow(a1[0],p-2,p))%p
                if not union(a1[1],a2[1],c): conf+=1
                nr+=1
            elif a1[1] is not None and a1[0] and (a2[1] is None or a2[0]==0):
                pass
        elif L[0]=='*':
            # b*(w - C) shape:  one factor a known constant
            X,Y=L[1],L[2]
            for A,B in ((X,Y),(Y,X)):
                c=None
                if A[0]=='c': c=A[1]%p
                elif A[0]=='v' and A[1] in BCONESET: c=val[A[1]]%p
                if c is None: continue
                if c==0: break
                if B[0]=='-' and B[2][0]=='c':
                    r=lin(B[1])
                    if r and r[1] is not None and r[0]:
                        forced[r[1]]=(B[2][1]%p)*pow(r[0],p-2,p)%p
                break
    return dict(Z=Z,find=find,forced=forced,nl=nl,nr=nr,conf=conf,val=val,lin=lin)
