#!/usr/bin/env python3
"""AUDIT T20 -- reconcile the 26: my 3,707 atoms of shape (w - (P*u)) with P in the p-class,
against L's 3,681 residual-atom census (which L's incidence criterion and hand-off argument use)."""
import os,sys,pickle,collections,re
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from circ2 import vars_of
import engine as E
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
v2a=collections.defaultdict(list)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
par={}
def find(x):
    par.setdefault(x,x)
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
def uni(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: par[ra]=rb
cp=re.compile(r'^\(x(\d+)-x(\d+)\)$')
for a in names:
    m=cp.match(a.replace(' ',''))
    if m: uni(int(m.group(1)),int(m.group(2)))
PCLASS={x for x in par if find(x)==find(26064)}
prod=re.compile(r'^\(x(\d+)-\(x(\d+)\*x(\d+)\)\)$')
mine={}   # u -> (atom index, h, P)
dupu=collections.Counter()
for i,a in enumerate(names):
    m=prod.match(a.replace(' ',''))
    if not m: continue
    h,Pv,u=int(m.group(1)),int(m.group(2)),int(m.group(3))
    if Pv in PCLASS:
        mine[u]=(i,h,Pv); dupu[u]+=1
print('atoms (w - (P*u)) with P in p-class : %d'%sum(dupu.values()))
print('distinct cofactors u among them     : %d'%len(mine))
print('cofactors appearing more than once  : %d'%sum(1 for u in dupu if dupu[u]>1))
H=pickle.load(open(os.path.join(LAB,'agentL_work','handles.pkl'),'rb'))
Lu=set(H['handle'])
mu=set(mine)
print('\nL\'s census: %d      mine: %d'%(len(Lu),len(mu)))
print('  mine \\ L : %d'%len(mu-Lu))
print('  L \\ mine : %d'%len(Lu-mu))
extra=sorted(mu-Lu)
print('\n== the %d atoms L does not count =='%len(extra))
for u in extra:
    i,h,Pv=mine[u]
    others=[j for j in v2a[h] if j!=i]
    g=[names[j] for j in others]
    print('  u=x%-6d free=%-5s  h=x%-6d P=x%-6d  guards:%d  %s'%(
        u,E.definer[u] is None,h,Pv,len(others),(g[0][:52] if g else '(NO GUARD)')))
missing=sorted(Lu-mu)
if missing:
    print('\n== the %d L counts that I do not =='%len(missing))
    for u in missing[:30]:
        ai=v2a[u]
        print('  u=x%-6d in %d atoms: %s'%(u,len(ai),[names[j][:46] for j in ai][:2]))
