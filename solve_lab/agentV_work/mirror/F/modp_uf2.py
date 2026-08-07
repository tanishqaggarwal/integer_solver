#!/usr/bin/env python3
"""Mod-p rigidity, version 2:
 (0) compute Z = wires provably == 0 mod p for EVERY assignment (closure over the definition DAG);
 (1) weighted union-find over wires using scalar defs, additive defs with a Z-summand, and the
     rigid residual atoms;
 (2) report which class each conditional-pin wire lands in."""
import sys,os,re,json,pickle,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663

# ---------- (0) Z-closure ----------
defrhs={}
for a in E.order:
    k,v,rhs=E.cls[a]; defrhs[v]=rhs
Z=set()
def isZ(n):
    o=n[0]
    if o=='c': return n[1]%p==0
    if o=='v': return n[1] in Z
    if o=='neg': return isZ(n[1])
    if o=='*': return isZ(n[1]) or isZ(n[2])
    return isZ(n[1]) and isZ(n[2])          # + and -
changed=True
while changed:
    changed=False
    for v,rhs in defrhs.items():
        if v in Z: continue
        if isZ(rhs): Z.add(v); changed=True
print('wires provably == 0 (mod p) for every assignment:',len(Z))

# ---------- (1) weighted union-find ----------
par=list(range(NV)); fac=[1]*NV
def find(v):
    if par[v]==v: return v,1
    r,f=find(par[v]); par[v]=r; fac[v]=fac[v]*f%p; return r,fac[v]
def union(a,b,c):
    ra,fa=find(a); rb,fb=find(b)
    if ra==rb: return (fa*pow(fb,p-2,p))%p==c%p
    par[ra]=rb; fac[ra]=c%p*fb%p*pow(fa,p-2,p)%p
    return True

def lin(n):
    """return (coef, var) if n reduces mod p to coef*x_var, else None."""
    o=n[0]
    if o=='v': return (0,None) if n[1] in Z else (1,n[1])
    if o=='c': return (0,None) if n[1]%p==0 else None
    if o=='neg':
        r=lin(n[1]); return None if r is None else ((-r[0])%p, r[1])
    if o=='*':
        a=lin(n[1]); b=lin(n[2])
        if a is None or b is None: return None
        if a[1] is None or a[0]==0: return (0,None)
        if b[1] is None or b[0]==0: return (0,None)
        return None                       # genuine product of two live wires
    a=lin(n[1]); b=lin(n[2])
    if a is None or b is None: return None
    sg=1 if o=='+' else -1
    if a[1] is None or a[0]==0: return ((sg*b[0])%p, b[1])
    if b[1] is None or b[0]==0: return (a[0]%p, a[1])
    if a[1]==b[1]: 
        c=(a[0]+sg*b[0])%p
        return (c, a[1]) if c else (0,None)
    return None

nl=0; conf=0
for v,rhs in defrhs.items():
    if v in Z: continue
    r=lin(rhs)
    if r is None or r[1] is None: continue
    if not union(v,r[1],r[0]): conf+=1
    nl+=1
print('rigid links from definitions:',nl,'conflicts',conf)

pats=[re.compile(r'^\(\((.+)\)-\(\d+\*x(\d+)\)\)$'), re.compile(r'^\(\((.+)\)-x(\d+)\)$'),
      re.compile(r'^\(\((.+)\)\+x(\d+)\)$')]
used=0; conf2=0
for a in E.res:
    ast=E.atoms[a]
    # atom == 0 ; strip a trailing handle term that is in Z
    if ast[0] in ('-','+'):
        L,R=ast[1],ast[2]
        rl=lin(R)
        if rl is not None and (rl[1] is None or rl[0]==0):
            ll=lin(L)
            if ll is not None and ll[1] is not None:
                continue     # single wire forced 0 mod p; not a link
            # L is a difference of two live wires?
            if L[0] in ('-','+'):
                a1=lin(L[1]); a2=lin(L[2])
                if a1 and a2 and a1[1] is not None and a2[1] is not None:
                    sg=1 if L[0]=='+' else -1
                    # a1[0]*u + sg*a2[0]*v = 0  ->  u = -sg*a2/a1 * v
                    c=(-sg*a2[0]*pow(a1[0],p-2,p))%p
                    if not union(a1[1],a2[1],c): conf2+=1
                    used+=1
            elif L[0]=='*':
                pass
print('rigid links from residual atoms:',used,'conflicts',conf2)

# ---------- (2) pin wires vs coordinates ----------
pins=json.load(open(os.path.join(HERE,'pins.json')))
sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
A=set(sup['7715']); B=set(sup['34554'])
coord={'x1':12186,'y1':16742,'x2':14853,'y2':24908}
croot={k:find(v)[0] for k,v in coord.items()}
print('coordinate class roots:',croot)
res=collections.defaultdict(list)
for bit,lst in pins.items():
    tree='A' if int(bit) in A else 'B'
    for w,C in lst:
        r,f=find(w)
        which=tuple(sorted(k for k in croot if croot[k]==r))
        res[(tree,which)].append((int(bit),w,f,C))
for key,lst in sorted(res.items(), key=lambda kv:-len(kv[1])):
    tree,which=key
    vals=collections.defaultdict(list)
    for bit,w,f,C in lst:
        vals[C%p*pow(f,p-2,p)%p].append(bit)
    print('tree %s -> coordinate class %-14s : %3d pin wires, %3d distinct forced values'%(tree,which,len(lst),len(vals)))
