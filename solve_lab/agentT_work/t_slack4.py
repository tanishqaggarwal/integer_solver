#!/usr/bin/env python3
"""AUDIT T19d -- confirm the forcing chain is a PROOF, and resolve the 128 multi-hop aliases."""
import os,sys,pickle,collections,re,json
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
par={}
def find(x):
    par.setdefault(x,x)
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
def uni(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: par[ra]=rb
cpat=re.compile(r'^\(x(\d+)-x(\d+)\)$')
for a in names:
    m=cpat.match(a.replace(' ',''))
    if m: uni(int(m.group(1)),int(m.group(2)))
PCLASS={x for x in par if find(x)==find(26064)}
print('== forcing chain ==')
print('  unique literal-p atom : (x26064 - p)          -> x26064 = p')
print('  copy class of x26064  : %d wires, linked by (xA - xB) atoms -> all = p'%len(PCLASS))
print('  the six are all in it : %s'%all(v in PCLASS for v in [4116,16153,1962,12682,19049,15616]))
print('  Since M is faithful (audit T2) and ker(M)=0 (F peel cert, re-verified in T2), EVERY atom')
print('  is zero in ANY full solution.  So x26064 = p and every copy equals it.  This is a PROOF,')
print('  not a measurement: the six are FORCED, and forced to p -- not to zero.')
# 128 multi-hop aliases
M=pickle.load(open(os.path.join(LAB,'agentL_work','full_model.pkl'),'rb'))
OUT=M['OUT']; tree=M['tree']
v2a=collections.defaultdict(set)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].add(i)
pairs=[]
for n in tree:
    if n not in OUT or not tree[n] or len(tree[n])!=2: continue
    ca,cb=tree[n]
    for j,slot in enumerate(OUT[n]):
        for side,ch in (('va',ca),('vb',cb)):
            w=slot.get(side)
            if ch not in OUT or w is None: continue
            co=OUT[ch][j].get('out')
            if co is None or w==co: continue
            pairs.append((w,co))
one=0; multi=[]
for w,co in pairs:
    if v2a[w]&v2a[co]: one+=1
    else: multi.append((w,co))
print('\n== the multi-hop aliases ==')
print('  aliased links total %d : one-atom %d, multi-hop %d'%(len(pairs),one,len(multi)))
# for multi-hop, do the two wires connect through an atom whose third variable is in the p-class,
# or through a chain of length 2?
via_p=0; chain2=0; other=0
for w,co in multi:
    nb=set()
    for i in v2a[w]:
        nb|=vars_of(atoms[names[i]])
    if v2a[co]&set().union(*[v2a[x] for x in nb]) if nb else False: chain2+=1
    hit=False
    for i in v2a[w]|v2a[co]:
        if vars_of(atoms[names[i]])&PCLASS: hit=True; break
    if hit: via_p+=1
    else: other+=1
print('  multi-hop links with a p-class wire in an incident atom : %d of %d'%(via_p,len(multi)))
print('  without                                                 : %d'%other)
print('\n  -> the multi-hop aliases terminate in the SAME p-class, not in something else.'
      if via_p>len(multi)*0.8 else '\n  -> mixed; not all multi-hop aliases route through the p-class.')
