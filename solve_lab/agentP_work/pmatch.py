#!/usr/bin/env python3
"""Agent P: orient atoms into an SLP by bipartite matching + acyclicity check."""
import pickle,sys,json
from collections import defaultdict,Counter,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
rows,AP=D['rows'],D['AP']
NV=38748

def cands(ap):
    """vars appearing ONLY in a lone degree-1 monomial."""
    cnt=Counter()
    for m in ap:
        for x in set(m): cnt[x]+=1
        for x in m:
            pass
    # count monomials each var occurs in, and whether it has a lone linear monomial
    occ=defaultdict(int); lin=set(); nonlin=set()
    for m in ap:
        s=set(m)
        for x in s: occ[x]+=1
        if len(m)==1: lin.add(m[0])
        else:
            for x in s: nonlin.add(x)
    return [x for x in lin if occ[x]==1 and x not in nonlin]

C=[cands(ap) for ap in AP]
print("cand-count hist:",sorted(Counter(len(c) for c in C).items()))

var2at=defaultdict(list)
for i,c in enumerate(C):
    for x in c: var2at[x].append(i)

# Hopcroft-Karp style simple Kuhn matching: match atoms -> vars (maximize matched vars)
matchA=[-1]*len(AP)   # atom -> var
matchV={}             # var -> atom
def try_k(a,seen):
    for x in C[a]:
        if x in seen: continue
        seen.add(x)
        if x not in matchV or try_k(matchV[x],seen):
            matchV[x]=a; matchA[a]=x; return True
    return False
sys.setrecursionlimit(300000)
order=sorted(range(len(AP)), key=lambda a: len(C[a]))
n=0
for a in order:
    if C[a] and try_k(a,set()): n+=1
print("matched atoms:",n,"of",len(AP),"; defined vars:",len(matchV))
free=[x for x in range(NV) if x not in matchV]
print("undefined (free) vars:",len(free))

# acyclicity check on the oriented DAG
indeg=[0]*len(AP)
deps=defaultdict(list)  # var -> atoms consuming it
defatom={x:a for x,a in matchV.items()}
for a,ap in enumerate(AP):
    if matchA[a]<0: continue
    ins=set()
    for m in ap: ins.update(m)
    ins.discard(matchA[a])
    for x in ins:
        if x in defatom: indeg[a]+=1; deps[x].append(a)
q=deque(a for a in range(len(AP)) if matchA[a]>=0 and indeg[a]==0)
topo=[]
while q:
    a=q.popleft(); topo.append(a)
    for b in deps[matchA[a]]:
        indeg[b]-=1
        if indeg[b]==0: q.append(b)
nm=sum(1 for a in range(len(AP)) if matchA[a]>=0)
print("topo sorted:",len(topo),"of",nm,"matched -> ACYCLIC" if len(topo)==nm else "-> HAS CYCLES")
pickle.dump({'matchA':matchA,'matchV':matchV,'topo':topo,'free':free,'C':C},open(W+'dag.pkl','wb'))
