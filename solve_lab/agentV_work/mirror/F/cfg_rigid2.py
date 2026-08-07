#!/usr/bin/env python3
"""Configuration-conditional mod-p rigidity, v2: AFFINE weighted union-find + general linear forms.

Fix a boolean configuration.  Wires whose free-input support lies inside the boolean free inputs have
values determined by that configuration alone, so selector products become linear with known constants.
Every relation derived below holds for EVERY integer assignment with that boolean configuration.
Representation:  value[x] = A[x] * value[root(x)] + B[x]   (mod p);  roots may additionally be pinned
to a constant in CONST.
"""
import sys,os,re,json,pickle,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
SUP=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
BOOLFREE=set()
for a in E.res:
    m=re.fullmatch(r'\(\(x(\d+)\*x(\d+)\)-x(\d+)\)',a) or re.fullmatch(r'\(x(\d+)\*\(1-x(\d+)\)\)',a) \
      or re.fullmatch(r'\(x(\d+)\*\(x(\d+)-1\)\)',a)
    if m and len(set(m.groups()))==1: BOOLFREE.add(int(m.group(1)))
BOOLFREE&=set(E.free)
BCONESET={v for v in range(NV) if set(SUP[str(v)])<=BOOLFREE}
defrhs={}
for a in E.order:
    k,v,rhs=E.cls[a]; defrhs[v]=rhs
inv=lambda a: pow(a%p,p-2,p)

def build(on_bits,verbose=False,rounds=4):
    v=[0]*NV
    for b in on_bits: v[b]=1
    E.run(v)
    val={w:v[w]%p for w in BCONESET}
    ZERO=set(w for w in BCONESET if val[w]==0)
    for rd in range(rounds):
        res=_pass(val,ZERO,verbose)
        newZ=set(ZERO)
        rel=res['rel']; CONST=res['CONST']
        for w in range(NV):
            r,a,b=rel(w)
            if r in CONST and (a*CONST[r]+b)%p==0: newZ.add(w)
        if verbose: print('   round %d: |ZERO| %d -> %d, links %d/%d, conflicts %d'%(rd,len(ZERO),len(newZ),res['nl'],res['nr'],res['conflicts']),flush=True)
        if newZ==ZERO: break
        ZERO=newZ
    res['ZERO']=ZERO
    return res

def _pass(val,ZERO,verbose=False):
    par=list(range(NV)); A=[1]*NV; B=[0]*NV
    CONST={}
    def find(x):
        path=[]
        while par[x]!=x: path.append(x); x=par[x]
        a,b=1,0
        for y in reversed(path):
            a,b=A[y]*a%p,(A[y]*b+B[y])%p
            par[y]=x; A[y]=a; B[y]=b
        return x
    def get(x):
        r=find(x); return r,A[x] if par[x]!=x else 1, B[x] if par[x]!=x else 0
    def rel(x):
        r=find(x)
        return (r,1,0) if r==x else (r,A[x],B[x])
    conflicts=[0]
    def union(u,c,vv,d):
        """assert value[u] = c*value[vv] + d  (mod p)"""
        ru,au,bu=rel(u); rv,av,bv=rel(vv)
        c%=p; d%=p
        # au*val(ru)+bu = c*(av*val(rv)+bv)+d
        if ru==rv:
            # (au - c*av) val(ru) = c*bv + d - bu
            k=(au-c*av)%p; m=(c*bv+d-bu)%p
            if k==0:
                if m!=0: conflicts[0]+=1
                return
            newc=m*inv(k)%p
            if ru in CONST and CONST[ru]!=newc: conflicts[0]+=1
            CONST[ru]=newc; return
        na=c*av%p*inv(au)%p; nb=(c*bv+d-bu)%p*inv(au)%p
        par[ru]=rv; A[ru]=na; B[ru]=nb
        if ru in CONST:
            cu=CONST.pop(ru)          # val(ru)=cu  ->  val(rv) determined
            if na: 
                cv=(cu-nb)%p*inv(na)%p
                if rv in CONST and CONST[rv]!=cv: conflicts[0]+=1
                CONST[rv]=cv
    def lf(n):
        """linear form: (dict var->coef, const) mod p, or None if genuinely nonlinear."""
        o=n[0]
        if o=='v':
            w=n[1]
            if w in ZERO: return ({},0)
            if w in BCONESET: return ({},val[w])
            return ({w:1},0)
        if o=='c': return ({},n[1]%p)
        if o=='neg':
            r=lf(n[1]); return None if r is None else ({k:(-x)%p for k,x in r[0].items()},(-r[1])%p)
        if o in '+-':
            a=lf(n[1]); b=lf(n[2])
            if a is None or b is None: return None
            s=1 if o=='+' else -1
            d=dict(a[0])
            for k,x in b[0].items():
                d[k]=(d.get(k,0)+s*x)%p
                if d[k]==0: del d[k]
            return (d,(a[1]+s*b[1])%p)
        a=lf(n[1]); b=lf(n[2])
        if a is None or b is None: return None
        if not a[0]: return ({k:a[1]*x%p for k,x in b[0].items()} if a[1] else {}, a[1]*b[1]%p)
        if not b[0]: return ({k:b[1]*x%p for k,x in a[0].items()} if b[1] else {}, b[1]*a[1]%p)
        return None
    def apply(form,extra_const=0):
        """assert  form == extra_const  (mod p)"""
        if form is None: return 0
        d,c=form; c=(c-extra_const)%p
        ks=list(d)
        if len(ks)==0:
            if c: conflicts[0]+=1
            return 0
        if len(ks)==1:
            u=ks[0]
            # d[u]*u + c = 0  -> u = -c/d[u]
            ru,au,bu=rel(u); target=(-c)%p*inv(d[u])%p
            # au*val(ru)+bu = target
            if au==0: return 0
            cv=(target-bu)%p*inv(au)%p
            if ru in CONST and CONST[ru]!=cv: conflicts[0]+=1
            CONST[ru]=cv; return 1
        if len(ks)==2:
            u,w=ks
            union(u,(-d[w])%p*inv(d[u])%p,w,(-c)%p*inv(d[u])%p); return 1
        return 0
    # close ZERO through the definition DAG (products with a zero factor, sums of zeros)
    def isZ(n):
        o=n[0]
        if o=='c': return n[1]%p==0
        if o=='v': return n[1] in ZERO or (n[1] in BCONESET and val[n[1]]==0)
        if o=='neg': return isZ(n[1])
        if o=='*': return isZ(n[1]) or isZ(n[2])
        return isZ(n[1]) and isZ(n[2])
    ch=True
    while ch:
        ch=False
        for w,rhs in defrhs.items():
            if w in ZERO: continue
            if isZ(rhs): ZERO.add(w); ch=True
    nl=0
    # definitions
    for w,rhs in defrhs.items():
        f=lf(rhs)
        if f is None: continue
        d,c=dict(f[0]),f[1]
        d[w]=(d.get(w,0)-1)%p
        if d[w]==0: del d[w]
        nl+=apply((d,c))
    # residual atoms (each must be 0)
    nr=0
    for a in E.res:
        f=lf(E.atoms[a])
        if f is None: continue
        nr+=apply(f)
    return dict(find=find,rel=rel,CONST=CONST,nl=nl,nr=nr,conflicts=conflicts[0],val=val,lf=lf,
                value_of=lambda x:(lambda r: CONST.get(r[0]) if r[0] in CONST else None)(rel(x)),
                rel_of=rel)
