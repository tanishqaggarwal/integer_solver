#!/usr/bin/env python3
"""AUDIT T19 -- the last gate: is anything in the instance forcing the six shared slack factors
to zero?   x_4116, x_16153, x_1962, x_12682, x_19049, x_15616
Answered from agent F's 39,033-atom decomposition, certified faithful in audit T2, so this is an
independent route from L's calibrated model."""
import os,sys,pickle,collections,re,json
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from circ2 import vars_of
import checker as CK, engine as E
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
v2a=collections.defaultdict(list)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
a2e=collections.defaultdict(set)
for e,row in enumerate(eqrows):
    for k,a in row: a2e[idx[a]].add(e)
codes,varsets=CK.load_equations()
v2e=collections.defaultdict(set)
for e,vs in enumerate(varsets):
    for u in vs: v2e[u].add(e)
SIX=[4116,16153,1962,12682,19049,15616]
print('='*90)
for v in SIX:
    A=v2a[v]
    print('\nx_%d :  free=%s   in %d atoms, %d equations'%(v,E.definer[v] is None,len(A),len(v2e[v])))
    unary=[]; boolp=[]; other=[]
    for i in A:
        s=names[i].replace(' ','')
        vs=vars_of(atoms[names[i]])
        if vs=={v}:
            unary.append(s)
            if re.fullmatch(r'\(?x%d\*\(x%d-1\)\)?'%(v,v),s) or re.fullmatch(r'\(?\(x%d\*x%d\)-x%d\)?'%(v,v,v),s): boolp.append(s)
        else: other.append(s)
    print('   atoms involving ONLY x_%d (a unary pin): %d  %s'%(v,len(unary),unary[:4]))
    print('   of those, boolean pins x*(x-1): %d'%len(boolp))
    print('   sample of the other atoms:')
    for s in other[:6]: print('      %s'%s[:110])
    # how does v occur -- always as a factor in a product?
    prod=sum(1 for i in A if re.search(r'\*x%d(?![0-9])'%v, names[i].replace(' ','')) or re.search(r'x%d(?![0-9])\*'%v, names[i].replace(' ','')))
    print('   occurrences where x_%d is a FACTOR in a product: %d of %d'%(v,prod,len(A)))
print('\n'+'='*90)
print('SUMMARY: any of the six carrying a unary pin?')
for v in SIX:
    un=[i for i in v2a[v] if vars_of(atoms[names[i]])=={v}]
    print('   x_%-6d unary-pin atoms: %d'%(v,len(un)))
