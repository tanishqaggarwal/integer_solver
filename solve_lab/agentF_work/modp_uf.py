#!/usr/bin/env python3
"""Mod-p rigid-link analysis: weighted union-find over wires, from (a) scalar definitions and
(b) residual atoms whose mod-p reduction is  alpha*u = beta*v  with alpha,beta units."""
import sys,os,re,json,pickle,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663

par=list(range(NV)); fac=[1]*NV   # value[v] = fac[v] * value[find(v)]  (mod p)
def find(v):
    if par[v]==v: return v,1
    r,f=find(par[v]); par[v]=r; fac[v]=fac[v]*f%p; return r,fac[v]
def union(a,b,c):
    """assert  value[a] = c * value[b]  (mod p)"""
    ra,fa=find(a); rb,fb=find(b)
    if ra==rb: return (fa*pow(fb,p-2,p))%p==c%p
    # fa*val(ra) = c*fb*val(rb)  ->  val(ra) = c*fb/fa * val(rb)
    par[ra]=rb; fac[ra]=c%p*fb%p*pow(fa,p-2,p)%p
    return True

# (a) scalar definitions  x_a := x_b   or   x_a := C*x_b   or   x_a := x_b*C
nlink=0; bad=0
for a in E.order:
    k,v,rhs=E.cls[a]
    s=node_str(rhs)
    m=re.fullmatch(r'x(\d+)',s)
    if m: nlink+= union(v,int(m.group(1)),1); continue
    m=re.fullmatch(r'\((\d+)\*x(\d+)\)',s) or re.fullmatch(r'\(x(\d+)\*(\d+)\)',s)
    if m:
        g1,g2=m.group(1),m.group(2)
        if g1.isdigit() and not g2.isdigit(): pass
        try:
            C=int(g1); w=int(g2)
        except ValueError:
            C=int(g2); w=int(g1)
        if s.startswith('(x'): C=int(g2); w=int(g1)
        nlink+= union(v,w,C); continue
print('scalar-definition links:',nlink)

# (b) residual atoms giving alpha*u = beta*v mod p (handles are 0 mod p)
pats=[ (re.compile(r'^\(\(x(\d+)-x(\d+)\)-\((\d+)\*x(\d+)\)\)$'), 'uv_Mh'),
       (re.compile(r'^\(\(x(\d+)-x(\d+)\)-x(\d+)\)$'),           'uv_h'),
       (re.compile(r'^\(\(x(\d+)-x(\d+)\)\+x(\d+)\)$'),          'uv_ph'),
       (re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-x(\d+)\)$'),'Cuv_h'),
       (re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-\((\d+)\*x(\d+)\)\)$'),'Cuv_Mh') ]
cnt=collections.Counter(); conflicts=0; used=0
for a in E.res:
    for pat,tag in pats:
        m=pat.match(a)
        if not m: continue
        g=m.groups()
        if tag in ('uv_Mh','uv_h','uv_ph'): u,v=int(g[0]),int(g[1])
        else: u,v=int(g[1]),int(g[2])
        if not union(u,v,1): conflicts+=1
        cnt[tag]+=1; used+=1
        break
print('rigid mod-p equality atoms used:',used,dict(cnt),'conflicts:',conflicts)

# now: classes of the pin wires
pins=json.load(open(os.path.join(HERE,'pins.json')))
sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
A=set(sup['7715']); B=set(sup['34554'])
coord={'x1':12186,'y1':16742,'x2':14853,'y2':24908,'x3':22162,'y3':30213}
croot={k:find(v)[0] for k,v in coord.items()}
print('coordinate roots:',{k:croot[k] for k in croot})
byclass=collections.defaultdict(list)
for bit,lst in pins.items():
    tree='A' if int(bit) in A else ('B' if int(bit) in B else '?')
    for w,C in lst:
        r,f=find(w)
        which=[k for k in croot if croot[k]==r]
        byclass[(tree,r,tuple(which))].append((int(bit),w,f,C))
print()
for (tree,r,which),lst in sorted(byclass.items(), key=lambda kv:-len(kv[1])):
    print('tree %s  class root x%-6d  contains coordinate(s) %s   ->  %d pin wires'%(tree,r,which,len(lst)))
# rigidity check: within a class, all pins force  coord = f^-1 * C ; distinct bits must agree
print()
for (tree,r,which),lst in byclass.items():
    if not which: continue
    vals={}
    for bit,w,f,C in lst:
        val=C%p*pow(f,p-2,p)%p
        vals.setdefault(val,[]).append(bit)
    print('class %s/%s : %d bits, %d DISTINCT forced values for %s  => at most one bit of this class can be ON'%(
        tree,which,len(lst),len(vals),which))
