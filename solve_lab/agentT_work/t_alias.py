#!/usr/bin/env python3
"""AUDIT T18 -- item 2: was agent K's withdrawal correct, or did K's null measure the aliasing?
Q reports no slot output feeds another slot directly (0/383), with an additive alias layer
  x_24468 = x_13682 + 12354891*x_34243
in between.  If so, a search for a DIRECT composition would find nothing even though the
composition exists through the alias.  Third independent check, from L's calibrated model."""
import os,sys,pickle,collections,re
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F)
from circ2 import vars_of
M=pickle.load(open(os.path.join(LAB,'agentL_work','full_model.pkl'),'rb'))
OUT=M['OUT']; tree=M['tree']
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
v2a=collections.defaultdict(list)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
direct=0; indirect=0; missing=0; pairs=[]
for n in tree:
    if not tree[n] or len(tree[n])!=2: continue
    ca,cb=tree[n]
    if n not in OUT or not tree[n]: continue
    for j,slot in enumerate(OUT[n]):
        for side,ch in (('va',ca),('vb',cb)):
            w=slot.get(side)
            if ch not in OUT or w is None: missing+=1; continue
            co=OUT[ch][j].get('out')
            if co is None: missing+=1; continue
            if w==co: direct+=1
            else:
                indirect+=1
                if len(pairs)<400: pairs.append((w,co))
tot=direct+indirect
print('parent-slot / child-output links examined: %d  (unresolvable: %d)'%(tot,missing))
print('   child out IS the parent slot wire (DIRECT) : %d'%direct)
print('   child out is a DIFFERENT wire (aliased)    : %d'%indirect)
print('   -> direct rate %.1f%%'%(100.0*direct/max(tot,1)))
# for the aliased ones, is there an atom containing BOTH wires?
lin=0; shapes=collections.Counter()
for w,co in pairs:
    common=set(v2a[w])&set(v2a[co])
    if common:
        lin+=1
        s=names[sorted(common)[0]]
        shapes[re.sub(r'x\d+','x',s.replace(' ',''))]+=1
print('\nof %d sampled aliased links, %d have an atom containing BOTH wires'%(len(pairs),lin))
print('commonest alias atom shapes:')
for s,n in shapes.most_common(6): print('   %-40s %d'%(s,n))
