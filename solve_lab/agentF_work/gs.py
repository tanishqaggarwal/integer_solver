#!/usr/bin/env python3
import sys,os,json,time,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from repair import *
def run(seed_assign, maxsweep=40, log=True):
    v=[0]*NV
    for k,x in seed_assign.items(): v[k]=x
    hist=[]
    for sweep in range(maxsweep):
        r,bad=full(v)
        nz=nz_atoms(r)
        if log: print('sweep',sweep,'nonzero',len(nz),'score',39033-len(bad),flush=True)
        hist.append((len(nz),39033-len(bad)))
        if not nz: break
        prog=False
        for i in list(nz):
            r=E.run(v)
            if r[i]==0: continue
            res=try_fix(v,i)
            if res: prog=True
        if not prog:
            if log: print('  no progress'); 
            break
    r,bad=full(v)
    return v,39033-len(bad),nz_atoms(r)
if __name__=='__main__':
    K1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    K2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    t0=time.time()
    v,s,nz=run({22106:1,5090:1,22162:K1,30213:K2,24468:K1,18956:K2})
    print('FINAL score',s,'nonzero',len(nz),'t',time.time()-t0)
    for i in nz: print('   ',E.res[i][:130])
    json.dump({'x_%d'%i:v[i] for i in range(NV) if v[i]},open(os.path.join(HERE,'gs1.json'),'w'))
