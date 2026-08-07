import sys, os, json, collections, pickle
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
# Z closure: vars that are 0 mod p for every assignment
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
print('|Z| =',len(Z))
# handle free vars: free vars u such that some Z-var is p*u style; find free vars appearing ONLY under Z
freevars=set(range(NV))-set(defrhs)
# which free vars appear (transitively) only inside Z-vars
appear=collections.defaultdict(set)  # free var -> set of atoms it influences
# build reverse: atom -> free vars via DAG
memo={}
def freeof(v):
    if v in memo: return memo[v]
    if v not in defrhs: memo[v]={v}; return memo[v]
    s=set()
    for u in vars_of(defrhs[v]): s|=freeof(u)
    memo[v]=s; return s
atomfree={}
for a in E.res:
    s=set()
    for u in vars_of(E.atoms[a]): s|=freeof(u)
    atomfree[a]=s
    for u in s: appear[u].add(a)
# classify free vars
hz=[u for u in freevars if u in Z or all(False for _ in ())]
zfree=[u for u in freevars if u in Z]
print('free vars',len(freevars),' of which in Z (handle roots):',len(zfree))
# handle vars: free vars whose every path multiplies by p -> they are in Z? no, u itself isn't in Z.
# instead: find free vars u such that for each atom a containing u, u enters only via a Z-var
cnt=collections.Counter(len(appear[u]) for u in freevars)
print('free var -> #atoms hist',sorted(cnt.items())[:12])
# Now: identify handle free vars = those that occur in exactly 1 atom AND enter through a p-multiplied var
Zfreeof={}
def zfreeof(v):
    """free vars reachable only through Z-multiplied subterms"""
    pass
# simpler: a free var u is a HANDLE if setting u changes the atom value by a multiple of p
import random
handle=set()
nonh=set()
prog=[]
for a in E.order:
    c=E.cls[a]; prog.append(c)
v0=[random.randrange(1,1000) for _ in range(NV)]
def evalatoms(v):
    return E.run(list(v))
r0=evalatoms(v0)
print('probing...')
for u in sorted(freevars):
    if len(appear[u])==0: continue
    v1=list(v0); v1[u]+=1
    r1=evalatoms(v1)
    d=[(x-y) for x,y in zip(r1,r0)]
    nzi=[i for i,x in enumerate(d) if x]
    if nzi and all(x%p==0 for x in d):
        handle.add(u)
    else: nonh.add(u)
print('handle free vars (every atom-delta divisible by p):',len(handle))
print('non-handle free vars:',len(nonh))
cnt2=collections.Counter(len(appear[u]) for u in handle)
print('handle -> #atoms hist',sorted(cnt2.items()))
pickle.dump({'Z':Z,'handle':handle,'nonh':nonh,'appear':{u:sorted(appear[u]) for u in freevars}},open('handles.pkl','wb'))
