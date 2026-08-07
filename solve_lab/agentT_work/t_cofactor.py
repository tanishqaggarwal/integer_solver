#!/usr/bin/env python3
"""AUDIT T16 -- the premise under L's incidence map, and the h-vs-u distinction.

L's criterion: 'every residual atom has exactly ONE free cofactor u and u occurs nowhere else,
so equation e contains atom a  <=>  u_a in vars(e)'.  That criterion produced 'of 3,681 atoms
exactly 15 are incident to the 25', which is the filter M is enumerating against.
Test the premise across all 3,681, in agent F's decomposition (certified faithful in audit T2).

Also: the coordinator reports M pricing a lattice target on x642 and x28730 'two of the twelve
cofactor variables'.  They are NOT in L's twelve.  Check what they actually are."""
import os,sys,json,pickle,collections
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from circ2 import vars_of
import checker as CK, engine as E
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
av=[frozenset(vars_of(atoms[a])) for a in names]
v2a=collections.defaultdict(set)
for i,vs in enumerate(av):
    for u in vs: v2a[u].add(i)
a2e=collections.defaultdict(set)
for e,row in enumerate(eqrows):
    for k,a in row: a2e[idx[a]].add(e)
codes,varsets=CK.load_equations()
v2e=collections.defaultdict(set)
for e,vs in enumerate(varsets):
    for u in vs: v2e[u].add(e)
H=pickle.load(open(os.path.join(LAB,'agentL_work','handles.pkl'),'rb'))
U=sorted(set(H['handle']))
print("L's 'handle' list (the cofactors u): %d"%len(U))
free=sum(1 for u in U if E.definer[u] is None)
print('   free variables: %d of %d'%(free,len(U)))
cnt=collections.Counter(len(v2a[u]) for u in U)
print('\n== PREMISE 1: does each cofactor u occur in exactly ONE atom? (F\'s parse) ==')
for k in sorted(cnt): print('   u occurs in %d atom(s): %d cofactors'%(k,cnt[k]))
bad=[u for u in U if len(v2a[u])!=1]
print('   violations: %d'%len(bad))
print('\n== PREMISE 2: does  eqs(u) == eqs(atom_u)  exactly? (this is the criterion itself) ==')
ok=mismatch=noatom=0; examples=[]
for u in U:
    A=v2a[u]
    if len(A)!=1: noatom+=1; continue
    a=next(iter(A))
    if v2e[u]==a2e[a]: ok+=1
    else:
        mismatch+=1
        if len(examples)<5:
            examples.append((u,len(v2e[u]),len(a2e[a]),sorted(v2e[u]-a2e[a])[:4],sorted(a2e[a]-v2e[u])[:4]))
print('   exact match      : %d'%ok)
print('   MISMATCH         : %d'%mismatch)
print('   u not in exactly 1 atom (criterion undefined): %d'%noatom)
for u,nu,na,extra,miss in examples:
    print('     x%-6d eqs(u)=%d eqs(atom)=%d  u-only:%s atom-only:%s'%(u,nu,na,extra,miss))
print('\n== the h vs u distinction ==')
for v,lab in [(642,'x642'),(28730,'x28730'),(17325,'u of x642'),(9413,'u of x28730'),
              (31864,'h'),(10903,'u of x31864'),(29854,'h'),(1329,'u of x29854')]:
    print('   %-14s x%-6d free=%-5s occurs in %d atom(s), %d equation(s)'%(
        lab,v,E.definer[v] is None,len(v2a[v]),len(v2e[v])))
