#!/usr/bin/env python3
"""Produce a checkable certificate that ker(M)=0 over Z (and over any field of char > 80)."""
import numpy as np, scipy.sparse as sp, os, sys, time, json
HERE=os.path.dirname(os.path.abspath(__file__))
M=sp.load_npz(os.path.join(HERE,'M.npz')).tocsr()
n,m=M.shape
Mc=M.tocsc()
rowcols=[M.indices[M.indptr[i]:M.indptr[i+1]] for i in range(n)]
rowvals=[M.data[M.indptr[i]:M.indptr[i+1]] for i in range(n)]
colrows=[Mc.indices[Mc.indptr[j]:Mc.indptr[j+1]] for j in range(m)]
deg=np.array([len(r) for r in rowcols])
zero=np.zeros(m,dtype=bool)
order=[]
stack=[i for i in range(n) if deg[i]==1]
t0=time.time()
while stack:
    i=stack.pop()
    if deg[i]!=1: continue
    js=[c for c in rowcols[i] if not zero[c]]
    if len(js)!=1: continue
    j=int(js[0])
    if zero[j]: continue
    zero[j]=True; order.append((j,int(i)))
    for r in colrows[j]:
        deg[r]-=1
        if deg[r]==1: stack.append(int(r))
print('forced %d atoms; time %.1fs'%(len(order),time.time()-t0))
np.save(os.path.join(HERE,'peel_order.npy'),np.array(order,dtype=np.int64))
# ---- INDEPENDENT VERIFICATION (recomputed from scratch, no reuse of the state above) ----
M2=sp.load_npz(os.path.join(HERE,'M.npz')).tocsr()
ordr=np.load(os.path.join(HERE,'peel_order.npy'))
done=set(); ok=True
for k in range(len(ordr)):
    j,i=int(ordr[k,0]),int(ordr[k,1])
    cols=M2.indices[M2.indptr[i]:M2.indptr[i+1]]
    vals=M2.data[M2.indptr[i]:M2.indptr[i+1]]
    cj=None
    for c,v in zip(cols,vals):
        if c==j: cj=int(v)
        elif c not in done: ok=False; print('BAD step',k,'row',i,'atom',j,'col',c,'not yet zero'); break
    if cj is None or cj==0: ok=False; print('BAD step',k,'atom',j,'not in row',i); break
    if not ok: break
    done.add(j)
print('certificate verified:',ok,' atoms forced:',len(done),'of',m)
print('=> rank(M) =',m if ok and len(done)==m else 'unknown', ' dim ker(M) =',0 if ok and len(done)==m else 'unknown')
