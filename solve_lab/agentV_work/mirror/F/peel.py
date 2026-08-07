#!/usr/bin/env python3
"""Characteristic-free peeling: a row with a single surviving nonzero forces that atom to 0."""
import numpy as np, scipy.sparse as sp, collections, time, os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
M=sp.load_npz(os.path.join(HERE,'M.npz')).tocsr()
n,m=M.shape
Mc=M.tocsc()
rowcols=[M.indices[M.indptr[i]:M.indptr[i+1]].tolist() for i in range(n)]
colrows=[Mc.indices[Mc.indptr[j]:Mc.indptr[j+1]].tolist() for j in range(m)]
deg=np.array([len(r) for r in rowcols])
zero=np.zeros(m,dtype=bool)
t0=time.time()
stack=[i for i in range(n) if deg[i]==1]
nz=0
while stack:
    i=stack.pop()
    if deg[i]!=1: continue
    j=[c for c in rowcols[i] if not zero[c]]
    if len(j)!=1: continue
    j=j[0]
    if zero[j]: continue
    zero[j]=True; nz+=1
    for r in colrows[j]:
        deg[r]-=1
        if deg[r]==1: stack.append(r)
        elif deg[r]==0: pass
print('peeling forced %d of %d atoms to zero in %.1fs'%(nz,m,time.time()-t0))
surv=[j for j in range(m) if not zero[j]]
print('surviving atoms:',len(surv))
rows_left=[i for i in range(n) if deg[i]>0]
print('rows still carrying surviving atoms:',len(rows_left))
np.save(os.path.join(HERE,'peel_zero.npy'),zero)
