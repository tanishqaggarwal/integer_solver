#!/usr/bin/env python3
"""Chain propagation with knob freezing; allows non-affine knobs too."""
import sys,os,pickle,json,time,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from circ2 import vars_of
E=Engine(); NR=len(E.res)
sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
asup=[]
for a in E.res:
    s=set()
    for u in vars_of(E.atoms[a]): s|=set(sup[str(u)])
    asup.append(frozenset(s))
reach=collections.defaultdict(set)
for i,s in enumerate(asup):
    for f in s: reach[f].add(i)

def solve(v, maxit=200, verbose=True, frozen=None):
    if frozen is None: frozen=set()
    for it in range(maxit):
        r=E.run(v); nz=[i for i in range(NR) if r[i]]
        if verbose: print(' it',it,'nz',len(nz),flush=True)
        if not nz: return v,True
        prog=False
        for i in nz:
            r=E.run(v)
            if r[i]==0: continue
            base=r[i]
            cands=sorted(asup[i]-frozen, key=lambda f:(len(reach[f]),f))
            done=False
            for f in cands:
                old=v[f]
                v[f]=old+1; a1=E.run(v)[i]
                d=a1-base
                if d==0: v[f]=old; continue
                if base % d: v[f]=old; continue
                v[f]=old-base//d
                if E.run(v)[i]==0:
                    frozen.add(f); prog=True; done=True; break
                v[f]=old
            if not done and verbose: pass
        if not prog:
            if verbose: print('  stuck'); 
            return v,False
    return v,False

if __name__=='__main__':
    K1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    K2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    v=[0]*NV
    for k,x in {22106:1,5090:1,22162:K1,30213:K2,24468:K1,18956:K2}.items(): v[k]=x
    t0=time.time()
    v,ok=solve(v, frozen={22162,30213,24468,18956,22106,5090})
    r=E.run(v); bad=E.score(r)
    print('score',39033-len(bad),'nz',sum(1 for x in r if x),'ok',ok,'t',time.time()-t0)
    for i in range(NR):
        if r[i]: print('  NZ',E.res[i][:120],'=',str(r[i])[:30])
    json.dump({'x_%d'%i:v[i] for i in range(NV) if v[i]},open(os.path.join(HERE,'gs2.json'),'w'))
