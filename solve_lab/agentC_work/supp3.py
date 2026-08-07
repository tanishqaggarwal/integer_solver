import sys, collections, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
NV=L.NVARS
outs={}
for a,(c,t) in L.atom_out.items(): outs[t]=a
sccs=pickle.load(open('/home/user/integer_solver/solve_lab/agentC_work/sccs.pkl','rb'))
# break cycles: in each 2-SCC drop the definer of the second member -> becomes free
for c in sccs:
    if len(c)>1:
        v=max(c)          # make max free
        del outs[v]
free=[v for v in range(NV) if v not in outs]
print('free after cycle break',len(free))
# topo order
adj={v:([u for u in L.avars[outs[v]] if u!=v] if v in outs else []) for v in range(NV)}
indeg=collections.Counter()
users=collections.defaultdict(list)
for v in range(NV):
    for u in adj[v]:
        users[u].append(v); indeg[v]+=1
q=[v for v in range(NV) if indeg[v]==0]; topo=[]
while q:
    v=q.pop(); topo.append(v)
    for w in users[v]:
        indeg[w]-=1
        if indeg[w]==0: q.append(w)
print('topo covers',len(topo),'of',NV)
supp=[None]*NV
for v in topo:
    if v in outs:
        s=set()
        for u in adj[v]: s|=supp[u]
        supp[v]=frozenset(s)
    else: supp[v]=frozenset([v])
chk=[a for a in range(L.NA) if a not in L.atom_out]
csupp={}
for a in chk:
    s=set()
    for u in L.avars[a]: s|=supp[u]
    csupp[a]=frozenset(s)
hist=collections.Counter(len(s) for s in csupp.values())
print('supp hist:',sorted(hist.items())[:15],'max',max(len(s) for s in csupp.values()))
print('zero-support checks:',hist[0])
par=list(range(NV))
def f(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
for a,s in csupp.items():
    s=list(s)
    for i in range(1,len(s)):
        ra,rb=f(s[0]),f(s[i])
        if ra!=rb: par[ra]=rb
comp=collections.defaultdict(list)
for v in free: comp[f(v)].append(v)
sizes=sorted((len(c) for c in comp.values()),reverse=True)
print('components:',len(comp),'sizes',sizes[:15])
# checks per component
ccomp=collections.defaultdict(list)
for a,s in csupp.items():
    if s: ccomp[f(next(iter(s)))].append(a)
    else: ccomp['NONE'].append(a)
print('checks in biggest comp:',len(ccomp[max(comp,key=lambda k:len(comp[k]))]))
print('checks with no free support:',len(ccomp['NONE']))
pickle.dump({'outs':outs,'free':free,'supp':supp,'csupp':csupp,'topo':topo,'comp':dict(comp),'ccomp':dict(ccomp),'par':par},
            open('/home/user/integer_solver/solve_lab/agentC_work/supp3.pkl','wb'))
