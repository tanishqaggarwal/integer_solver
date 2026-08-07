import sys, os, json, collections, pickle
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
Z=set()
def isZ(n):
    o=n[0]
    if o=='c': return n[1]%p==0
    if o=='v': return n[1] in Z
    if o=='neg': return isZ(n[1])
    if o=='*': return isZ(n[1]) or isZ(n[2])
    return isZ(n[1]) and isZ(n[2])
ch=True
while ch:
    ch=False
    for v,r in defrhs.items():
        if v not in Z and isZ(r): Z.add(v); ch=True
memoP={}
def freePv(v):
    if v in memoP: return memoP[v]
    if v in Z: memoP[v]=frozenset(); return memoP[v]
    if v not in defrhs: memoP[v]=frozenset([v]); return memoP[v]
    memoP[v]=frozenset()  # cycle guard
    r=freeP(defrhs[v]); memoP[v]=r; return r
def freeP(n):
    o=n[0]
    if o=='c': return frozenset()
    if o=='v': return freePv(n[1])
    if o=='neg': return freeP(n[1])
    if o=='*':
        if isZ(n[1]) or isZ(n[2]): return frozenset()
        return freeP(n[1])|freeP(n[2])
    return freeP(n[1])|freeP(n[2])
memoF={}
def freeAll(v):
    if v in memoF: return memoF[v]
    if v not in defrhs: memoF[v]=frozenset([v]); return memoF[v]
    memoF[v]=frozenset()
    s=frozenset()
    for u in vars_of(defrhs[v]): s|=freeAll(u)
    memoF[v]=s; return s
sys.setrecursionlimit(100000)
atomP={}; atomA={}
for a in E.res:
    s=frozenset(); t=frozenset()
    for u in vars_of(E.atoms[a]): s|=freePv(u); t|=freeAll(u)
    atomP[a]=s; atomA[a]=t
freevars=set(range(NV))-set(defrhs)
appearP=collections.defaultdict(list); appearA=collections.defaultdict(list)
for a in E.res:
    for u in atomP[a]: appearP[u].append(a)
    for u in atomA[a]: appearA[u].append(a)
handle=[u for u in freevars if appearA[u] and not appearP[u]]
value=[u for u in freevars if appearP[u]]
print('free vars %d  |  handle-only %d  |  value %d  |  unused %d'%(len(freevars),len(handle),len(value),len(freevars)-len(handle)-len(value)))
print('handle -> #atoms hist',sorted(collections.Counter(len(appearA[u]) for u in handle).items()))
print('value  -> #atoms hist',sorted(collections.Counter(len(appearP[u]) for u in value).items()))
D=pickle.load(open('ortree2.pkl','rb')); tree=D['tree']
fl=set(v for v in tree if tree[v] is None and v not in defrhs)
grp=collections.defaultdict(list)
for u in value: grp[len(appearP[u])].append(u)
for k in sorted(grp):
    print(' value class %d atoms: %d vars; leaf-overlap %d'%(k,len(grp[k]),len(set(grp[k])&fl)))
pickle.dump({'Z':Z,'handle':handle,'value':value,'atomP':atomP,'appearP':{u:appearP[u] for u in value}},open('handles.pkl','wb'))
