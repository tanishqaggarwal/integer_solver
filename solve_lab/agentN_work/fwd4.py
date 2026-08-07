"""Extend the fwd2 orientation: let check atoms define currently-free vars when acyclic."""
import model, pickle, os, ast, re, time
from collections import defaultdict, Counter
HERE=os.path.dirname(os.path.abspath(__file__))
d=model.get(); atom_src=d['atom_src']; atom_vars=d['atom_vars']
cands=pickle.load(open(os.path.join(HERE,'defcands.pkl'),'rb'))
polys=pickle.load(open(os.path.join(HERE,'polys.pkl'),'rb'))
F=pickle.load(open(os.path.join(HERE,'fwd2.pkl'),'rb'))
definer=list(F['definer']); order=list(F['order']); free=set(F['free0'])
NV=38748; NA=len(atom_src)
used=set(a for a in definer if a>=0)

def compute(free_set, definer, order):
    fidx={v:i for i,v in enumerate(sorted(free_set))}
    sup=[0]*NV
    for v in free_set: sup[v]=1<<fidx[v]
    for v in order:
        a=definer[v]; s=0
        for u in atom_vars[a]:
            if u!=v: s|=sup[u]
        sup[v]=s
    return sup, fidx

sup,fidx=compute(free,definer,order)
print('start: free',len(free),'checks',NA-len(used))
for it in range(8):
    adopted=0
    # candidate check atoms that can define a currently free var
    for a in range(NA):
        if a in used: continue
        for v in cands[a]:
            if v not in free: continue
            others=[u for u in atom_vars[a] if u!=v]
            bit=1<<fidx[v]
            if any(sup[u]&bit for u in others): continue
            # adopt
            definer[v]=a; used.add(a); free.discard(v); order.append(v)
            adopted+=1
            break
    if adopted==0: break
    # re-topologically sort order by dependency (recompute)
    # rebuild order via Kahn over defined vars
    defs=[v for v in range(NV) if definer[v]>=0]
    indeg={}; children=defaultdict(list)
    for v in defs:
        deps=[u for u in atom_vars[definer[v]] if u!=v and definer[u]>=0]
        indeg[v]=len(deps)
        for u in deps: children[u].append(v)
    from collections import deque
    Q=deque(v for v in defs if indeg[v]==0); neworder=[]
    while Q:
        v=Q.popleft(); neworder.append(v)
        for w in children[v]:
            indeg[w]-=1
            if indeg[w]==0: Q.append(w)
    if len(neworder)!=len(defs):
        print('CYCLE detected at iter',it,len(neworder),len(defs)); break
    order=neworder
    sup,fidx=compute(free,definer,order)
    print('iter',it,'adopted',adopted,'free',len(free),'checks',NA-len(used))
pickle.dump({'definer':definer,'order':order,'free':sorted(free),'checks':sorted(set(range(NA))-used)},
            open(os.path.join(HERE,'fwd4.pkl'),'wb'))
print('FINAL free',len(free),'checks',NA-len(used))
