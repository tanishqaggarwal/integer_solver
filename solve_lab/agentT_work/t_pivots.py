#!/usr/bin/env python3
"""AUDIT T0: peel_cert.py verifies ker(M)=0 but NEVER checks pivot magnitudes -- it only
tests pivot != 0.  RESUME_F / FLEET.md claim 'all pivots are +-1 or +-2, none divisible by
any odd prime'.  Measure them."""
import os,sys,collections
import numpy as np, scipy.sparse as sp
F=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentF_work'))
M=sp.load_npz(os.path.join(F,'M.npz')).tocsr()
o=np.load(os.path.join(F,'peel_order.npy'))
piv=[]
for k in range(len(o)):
    j,i=int(o[k,0]),int(o[k,1])
    cols=list(M.indices[M.indptr[i]:M.indptr[i+1]]); vals=M.data[M.indptr[i]:M.indptr[i+1]]
    piv.append(int(vals[cols.index(j)]))
c=collections.Counter(piv)
print('peel steps:',len(piv))
print('distinct pivot values:',sorted(c.items()))
bad=[v for v in c if abs(v) not in (1,2)]
print("pivots outside {+-1,+-2}:",sorted(bad),' -> claim',('CONFIRMED' if not bad else 'FALSE'))
