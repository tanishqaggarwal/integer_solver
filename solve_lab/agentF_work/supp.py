#!/usr/bin/env python3
import sys,os,pickle,collections,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
free=set(E.free)
sup={}
for f in free: sup[f]=frozenset([f])
t0=time.time()
for a in E.order:
    c=E.cls[a]
    s=set()
    for u in vars_of(c[2]): s|=sup.get(u,frozenset())
    sup[c[1]]=frozenset(s)
print('supports built',time.time()-t0)
import statistics
sz=[len(sup[v]) for v in range(NV)]
print('support size: max',max(sz),'mean',sum(sz)/len(sz))
tgt={'x24468':24468,'x32989':32989,'x2300':2300,'x9274':9274,'x18956':18956,'x14257':14257}
for k,vv in tgt.items(): print(k,'support',len(sup[vv]))
pickle.dump({str(v):sorted(sup[v]) for v in range(NV)},open(os.path.join(HERE,'supp.pkl'),'wb'))
