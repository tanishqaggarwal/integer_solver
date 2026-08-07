#!/usr/bin/env python3
"""Agent P: recover global SLP order from intra-equation atom order, then orient."""
import pickle,sys,json
from collections import defaultdict,Counter,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
rows,AP=D['rows'],D['AP']
NA=len(AP); NV=38748

succ=defaultdict(set); indeg=Counter()
for r in rows:
    ids=[a for c,a in r['row']]
    for u,v in zip(ids,ids[1:]):
        if u!=v and v not in succ[u]:
            succ[u].add(v); indeg[v]+=1
print("order edges:",sum(len(s) for s in succ.values()))
q=deque(sorted(a for a in range(NA) if indeg[a]==0))
topo=[]
while q:
    a=q.popleft(); topo.append(a)
    for b in sorted(succ[a]):
        indeg[b]-=1
        if indeg[b]==0: q.append(b)
print("topo:",len(topo),"of",NA, "ACYCLIC" if len(topo)==NA else "CYCLIC")
if len(topo)!=NA:
    left=[a for a in range(NA) if a not in set(topo)]
    print("  unordered atoms:",len(left), left[:20])
    pickle.dump({'topo':topo,'succ':{k:list(v) for k,v in succ.items()}},open(W+'order.pkl','wb'))
    sys.exit()

# orient: walk in topo order, define the first unused candidate var
def cands(ap):
    occ=defaultdict(int); lin=set(); nonlin=set()
    for m in ap:
        for x in set(m): occ[x]+=1
        if len(m)==1: lin.add(m[0])
        else:
            for x in set(m): nonlin.add(x)
    return [x for x in lin if occ[x]==1 and x not in nonlin]
C=[cands(ap) for ap in AP]
defined={}
outof=[-1]*NA
constraints=[]
for a in topo:
    av=set()
    for m in AP[a]: av.update(m)
    undef=[x for x in C[a] if x not in defined]
    unknown=[x for x in av if x not in defined]
    if len(undef)==1:
        outof[a]=undef[0]; defined[undef[0]]=a
    elif len(undef)==0:
        constraints.append(a)
    else:
        # multiple undefined candidates -> pick the one that is NOT free-input-like:
        # heuristic: prefer the var with smallest number of later uses? just take max index of C order
        outof[a]=undef[-1]; defined[undef[-1]]=a
print("defined vars:",len(defined),"pure constraints:",len(constraints),"ambiguous handled")
free=[x for x in range(NV) if x not in defined]
print("free vars:",len(free))
pickle.dump({'topo':topo,'outof':outof,'defined':defined,'free':free,'C':C},open(W+'order.pkl','wb'))
