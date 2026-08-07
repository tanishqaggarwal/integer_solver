#!/usr/bin/env python3
"""Exact integer Jacobian of residual atoms wrt free inputs at a point; affineness test."""
import sys,os,pickle,json,time,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
E=Engine()
NR=len(E.res)

def build_jac(v0, knobs=None, check_affine=True):
    v=list(v0)
    F0=E.run(v)
    if knobs is None: knobs=E.free
    cols={}; nonaff=[]
    for f in knobs:
        old=v[f]
        v[f]=old+1; F1=E.run(v)
        col={i:F1[i]-F0[i] for i in range(NR) if F1[i]!=F0[i]}
        if check_affine and col:
            v[f]=old+2; F2=E.run(v)
            ok=all(F2[i]-F0[i]==2*col.get(i,0) for i in set(list(col)+[i for i in range(NR) if F2[i]!=F0[i]]))
            if not ok: nonaff.append(f); v[f]=old; continue
        v[f]=old
        if col: cols[f]=col
    return F0,cols,nonaff

if __name__=='__main__':
    K1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    K2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    v0=[0]*NV
    for k,x in {22106:1,5090:1,22162:K1,30213:K2,24468:K1,18956:K2}.items(): v0[k]=x
    t0=time.time()
    F0,cols,nonaff=build_jac(v0)
    print('time',time.time()-t0)
    print('knobs with nonzero column',len(cols),'non-affine knobs',len(nonaff))
    nz=[i for i in range(NR) if F0[i]]
    print('nonzero residual atoms',len(nz))
    print('total nnz',sum(len(c) for c in cols.values()))
    pickle.dump({'v0':v0,'F0':F0,'cols':cols,'nonaff':nonaff},open(os.path.join(HERE,'jac1.pkl'),'wb'))
