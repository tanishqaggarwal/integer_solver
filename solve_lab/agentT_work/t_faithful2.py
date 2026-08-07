#!/usr/bin/env python3
"""AUDIT T1b: faithfulness of M at MANY points, incl. random ones, vs checker.py's own evaluator."""
import sys,os,json,pickle,time,random
F=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentF_work'))
L=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
sys.path.insert(0,F); sys.path.insert(0,L)
import numpy as np, scipy.sparse as sp
from fwd import compile_node
import checker
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms); eqrows=d['eqrows']; idx={a:i for i,a in enumerate(names)}
NV=38748
src='r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']'
prog=compile(src,'<atoms>','exec')
t0=time.time(); eqs,_vs=checker.load_equations()
print('checker loaded %d eqs in %.0fs'%(len(eqs),time.time()-t0))

def model_bad(v):
    r=[0]*len(names); exec(prog,{'v':v,'r':r,'__builtins__':{}})
    bad=[]
    for e,row in enumerate(eqs if False else eqrows):
        s=0
        for k,a in row:
            val=r[idx[a]]
            if val: s+=k*val
        if s: bad.append(e)
    return bad,sum(1 for x in r if x)

def true_bad(v):
    return checker.evaluate_all(eqs,v)
def _old(v):
    env={'x_%d'%i:v[i] for i in range(NV)}
    bad=[]
    for e,f in enumerate(eqs):
        if f(env)!=0: bad.append(e)
    return bad

random.seed(7)
tests=[('zeros',[0]*NV)]
for name in ['new_instance_partial_39026','new_instance_partial_39024','new_instance_partial_39021','new_instance_partial_39013']:
    a=json.load(open(os.path.join(L,'best',name+'.json')))
    v=[0]*NV
    for k,val in a.items(): v[int(k[2:])]=int(val)
    tests.append((name,v))
for t in range(3):
    v=[random.randint(-5,5) for _ in range(NV)]
    tests.append(('random_small_%d'%t,v))
for t in range(2):
    v=[random.randint(-10**30,10**30) for _ in range(NV)]
    tests.append(('random_big_%d'%t,v))

for nm,v in tests:
    mb,nzA=model_bad(v)
    try:
        tb=true_bad(v)
        agree = (mb==tb)
        print('%-28s atoms_nz=%-6d model_bad=%-6d true_bad=%-6d AGREE=%s'%(nm,nzA,len(mb),len(tb),agree))
        if not agree:
            print('   only-model:',sorted(set(mb)-set(tb))[:10])
            print('   only-true :',sorted(set(tb)-set(mb))[:10])
    except Exception as ex:
        print('%-28s atoms_nz=%-6d model_bad=%-6d  TRUE EVAL FAILED: %s'%(nm,nzA,len(mb),type(ex).__name__))
