#!/usr/bin/env python3
"""AUDIT T19b -- full structure of the six shared slack factors: every atom, the definition,
and their values in the deliverable."""
import os,sys,pickle,collections,re,json
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from circ2 import vars_of
import engine as E, harness as H
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
v2a=collections.defaultdict(list)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
B=json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json')))
SIX=[4116,16153,1962,12682,19049,15616]
for v in SIX:
    print('\n'+'='*92)
    print('x_%d   definer=%s   deliverable value=%s'%(v,E.definer[v],B.get('x_%d'%v,'(absent=0)')))
    print('  ALL %d atoms containing it:'%len(v2a[v]))
    for i in v2a[v]:
        s=names[i].replace(' ','')
        tag=''
        if re.fullmatch(r'\(x\d+-\(x%d\*x\d+\)\)'%v,s) or re.fullmatch(r'\(x\d+-\(x\d+\*x%d\)\)'%v,s): tag='  [product w = v*z]'
        elif len(vars_of(atoms[names[i]]))==2 and '*' not in s: tag='  [copy/alias]'
        print('     %-72s%s'%(s[:72],tag))
