"""Greedy topological closure: repair every nonzero check by setting a free input."""
import sys, json, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
CHK=[a for a in range(L.NA) if a not in L.atom_out]
FREE=set(u for u in range(L.NVARS) if u not in outs)
TP={u:i for i,u in enumerate(topo)}
def lin_in(a,u,v):
    """(c, rest) with atom = c*v[u]+rest, or None if u nonlinear."""
    c=0; rest=0
    for m,cc in L.polys[a].items():
        k=m.count(u)
        if k>1: return None
        t=cc
        if k==0:
            for w in m: t*=v[w]
            rest+=t
        else:
            for w in m:
                if w!=u: t*=v[w]
            c+=t
    return c,rest
def evalat(a,v):
    s=0
    for m,c in L.polys[a].items():
        t=c
        for w in m: t*=v[w]
        s+=t
    return s
def closure(seeds, rounds=60, verbose=True):
    seeds=dict(seeds)
    hist=[]
    for r in range(rounds):
        v=[0]*L.NVARS
        for k,val in seeds.items(): v[k]=val
        forward(v)
        av=L.all_atom_values(v)
        nz=[a for a in CHK if av[a]!=0]
        f=L.failing_eqs(av)
        sc=L.NEQ-len(f)
        hist.append((sc,len(nz)))
        if verbose: print('  round %d: score=%d nz=%d nseeds=%d'%(r,sc,len(nz),len(seeds)))
        if not nz: return seeds,v,sc,nz
        prog=0
        # order: latest-first so early fixes are not undone
        nz.sort(key=lambda a: -max((TP[u] for u in L.avars[a]), default=0))
        for a in nz:
            if evalat(a,v)==0: continue
            cands=[]
            for u in L.avars[a]:
                if u not in FREE or u in seeds: continue
                r2=lin_in(a,u,v)
                if r2 is None: continue
                c,rest=r2
                if c==0 or rest%c: continue
                cands.append((TP[u],u,-rest//c))
            if not cands: continue
            cands.sort()
            _,u,val=cands[-1]
            seeds[u]=val; v[u]=val; forward(v); prog+=1
        if prog==0:
            if verbose: print('  STALL')
            break
    v=[0]*L.NVARS
    for k,val in seeds.items(): v[k]=val
    forward(v)
    av=L.all_atom_values(v)
    nz=[a for a in CHK if av[a]!=0]
    sc=L.NEQ-len(L.failing_eqs(av))
    return seeds,v,sc,nz
if __name__=='__main__':
    S={542:1, 91:1, 22162:K2, 30213:K1}
    t=time.time()
    seeds,v,sc,nz=closure(S)
    print('FINAL score',sc,'nz',len(nz),'time',time.time()-t)
    for a in nz[:25]:
        print('  a%-6d eqs=%-3d val=%s'%(a,len(L.atom2eq.get(a,{})),str(evalat(a,v))[:40]))
        print('     ',L.atom_src[a][:200])
    json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0}, open('/home/user/integer_solver/solve_lab/agentC_work/close_out.json','w'))
    json.dump({str(k):str(val) for k,val in seeds.items()}, open('/home/user/integer_solver/solve_lab/agentC_work/close_seeds.json','w'))
