#!/usr/bin/env python3
"""Build the 39,033 x 39,033 equation-atom incidence matrix M and save it."""
import sys,os,pickle,collections,time
import numpy as np, scipy.sparse as sp
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine
E=Engine()
names=list(E.atoms); idx={a:i for i,a in enumerate(names)}
rows=[];cols=[];vals=[]
for e,row in enumerate(E.eqrows):
    acc=collections.Counter()
    for k,a in row: acc[idx[a]]+=k          # merge repeated atoms within one equation
    for j,k in acc.items():
        if k: rows.append(e);cols.append(j);vals.append(k)
n=len(E.eqrows); m=len(names)
print('M is %d x %d with %d nonzeros'%(n,m,len(vals)))
M=sp.csr_matrix((np.array(vals,dtype=np.int64),(np.array(rows),np.array(cols))),shape=(n,m))
sp.save_npz(os.path.join(HERE,'M.npz'),M)
pickle.dump(names,open(os.path.join(HERE,'M_atomnames.pkl'),'wb'))
d_row=np.diff(M.indptr); Mc=M.tocsc(); d_col=np.diff(Mc.indptr)
print('row degree: min %d max %d mean %.2f ; zero rows %d'%(d_row.min(),d_row.max(),d_row.mean(),(d_row==0).sum()))
print('col degree: min %d max %d mean %.2f ; zero cols %d'%(d_col.min(),d_col.max(),d_col.mean(),(d_col==0).sum()))
print('rows of degree 1:',int((d_row==1).sum()),' cols of degree 1:',int((d_col==1).sum()))
print('max |coef|',int(np.abs(M.data).max()))
