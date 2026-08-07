#!/usr/bin/env python3
"""AUDIT T19c -- THE ANSWER.  The six 'shared slack factors' look like p-valued wires.
If so the question 'is anything forcing them to zero' is the wrong question: nothing does, and
nothing should -- they are forced to p, the slack term is p*u, and the alias is exact MOD P but
not over Z.  Establish it: (1) values, (2) the copy-equivalence class, (3) a literal pin to p."""
import os,sys,pickle,collections,re,json
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
B=json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json')))
SIX=[4116,16153,1962,12682,19049,15616]
print('== 1. deliverable values ==')
for v in SIX:
    val=B.get('x_%d'%v)
    print('   x_%-6d = %s'%(v,'EXACTLY p' if val is not None and int(val)==p else str(val)[:40]))
# copy-equivalence classes: atoms of shape (xA - xB)
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
print('\n== 2. copy-equivalence class of the six (atoms of shape (xA - xB)) ==')
cls=collections.defaultdict(list)
for x in list(par): cls[find(x)].append(x)
roots={find(v) for v in SIX if v in par}
print('   the six fall into %d copy class(es): %s'%(len(roots),sorted(roots)))
members=sorted(set().union(*[set(cls[r]) for r in roots]))
print('   total wires in that class: %d'%len(members))
inB=[m for m in members if B.get('x_%d'%m) is not None]
atp=[m for m in inB if int(B['x_%d'%m])==p]
print('   of those, set in the deliverable: %d ; equal to p: %d'%(len(inB),len(atp)))
print('\n== 3. is any member of the class pinned to the LITERAL p by an atom? ==')
lit=re.compile(r'^\(x(\d+)-(\d+)\)$')
hits=[]
for a in names:
    m=lit.match(a.replace(' ',''))
    if m and int(m.group(2))==p and int(m.group(1)) in set(members): hits.append(a)
print('   atoms of shape (x - p): %d  %s'%(len(hits),hits[:4]))
allp=[a for a in names if str(p) in a.replace(' ','')]
print('   atoms anywhere in the instance containing the literal p: %d'%len(allp))
for a in allp[:6]: print('      %s'%a[:100])
print('\n== 4. so what do the "slack" products look like? ==')
prod=re.compile(r'^\(x(\d+)-\(x(\d+)\*x(\d+)\)\)$')
n=0
for a in names:
    m=prod.match(a.replace(' ',''))
    if m and int(m.group(2)) in set(members):
        if n<5: print('      %-46s  ->  x%s = p * x%s'%(a[:46],m.group(1),m.group(3)))
        n+=1
print('   atoms of shape (w - (P*u)) with P in the p-class: %d'%n)
