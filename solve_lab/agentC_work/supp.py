import sys, collections, pickle, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
NV=L.NVARS
outs={}   # var -> atom defining it
for a,(c,t) in L.atom_out.items(): outs[t]=a
free=[v for v in range(NV) if v not in outs]
print('free',len(free))
# iterative support computation over topo order; cycles handled by fixpoint
supp=[None]*NV
for v in free: supp[v]=frozenset([v])
order=L.topo
defined_in_topo=set(order)
import sys as _s
_s.setrecursionlimit(100000)
# fixpoint
changed=True; rounds=0
cyc=[v for v in range(NV) if v in outs and v not in defined_in_topo]
print('cyclic vars',len(cyc))
for v in cyc: supp[v]=frozenset()
for v in order:
    a=outs[v]
    s=set()
    for u in L.avars[a]:
        if u==v: continue
        if supp[u] is None: supp[u]=frozenset()
        s|=supp[u]
    supp[v]=frozenset(s)
# now recompute cyclic by fixpoint
for r in range(30):
    ch=0
    for v in cyc:
        a=outs[v]; s=set()
        for u in L.avars[a]:
            if u==v: continue
            s|=supp[u]
        s=frozenset(s)
        if s!=supp[v]: supp[v]=s; ch+=1
    if ch==0: break
print('cyc rounds',r,'changed',ch)
chk=[a for a in range(L.NA) if a not in L.atom_out]
csupp={}
for a in chk:
    s=set()
    for u in L.avars[a]: s|=supp[u]
    csupp[a]=frozenset(s)
hist=collections.Counter(len(s) for s in csupp.values())
print('check free-support size hist (top 20):',sorted(hist.items())[:20])
print('max supp',max(len(s) for s in csupp.values()))
pickle.dump({'supp':supp,'csupp':csupp,'free':free,'outs':outs}, open('/home/user/integer_solver/solve_lab/agentC_work/supp.pkl','wb'))
# union-find over free inputs, joined by checks
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
print('n components over free inputs:',len(comp),'sizes top20:',sizes[:20])
