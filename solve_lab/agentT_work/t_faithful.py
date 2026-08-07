#!/usr/bin/env python3
"""AUDIT T1: is F's M a FAITHFUL linearisation?
Evaluate every atom at a GIVEN assignment (no forward re-derivation),
compute M*a, and compare the nonzero-row set to checker.py's failing set."""
import sys,os,json,pickle,time,collections
F=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentF_work')
F=os.path.abspath(F); sys.path.insert(0,F)
import numpy as np, scipy.sparse as sp
from fwd import compile_node

d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']
names=list(atoms)          # SAME order buildM.py used
NV=38748
assign=json.load(open(sys.argv[1]))
v=[0]*NV
for k,val in assign.items(): v[int(k[2:])]=int(val)

t0=time.time()
src='r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']'
prog=compile(src,'<atoms>','exec')
r=[0]*len(names)
exec(prog,{'v':v,'r':r,'__builtins__':{}})
print('atoms evaluated: %d  in %.1fs'%(len(names),time.time()-t0))
nz=[i for i,x in enumerate(r) if x]
print('nonzero atoms at this assignment:',len(nz))
print('  names:',[names[i] for i in nz][:30])

# M*a, exact python ints
M=sp.load_npz(os.path.join(F,'M.npz')).tocsr()
assert M.shape==(len(eqrows),len(names)), (M.shape,len(eqrows),len(names))
bad=[]
for e in range(M.shape[0]):
    s=0
    for j,k in zip(M.indices[M.indptr[e]:M.indptr[e+1]],M.data[M.indptr[e]:M.indptr[e+1]]):
        if r[j]: s+=int(k)*r[j]
    if s: bad.append(e)
print('rows with M*a != 0 :',len(bad))
print('  ',sorted(bad))

# also: directly from eqrows (independent of M) -- cross-check M against its source
idx={a:i for i,a in enumerate(names)}
bad2=[]
for e,row in enumerate(eqrows):
    s=0
    for k,a in row:
        val=r[idx[a]]
        if val: s+=k*val
    if s: bad2.append(e)
print('rows with eqrow-form != 0 :',len(bad2),' identical to M:',bad2==bad)
