"""Map M's 32 incident handles onto circuit roles, and ask which 4-subsets are realizable."""
import sys, json, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentF_work')
from fwd import Engine,NV
from circ2 import vars_of
from parse import node_str
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M32=[642,1627,1844,1956,2218,2892,4863,6480,7062,7945,9629,10861,11425,15422,16495,21279,
     21718,23538,23642,23754,23822,26732,28098,28730,29305,29854,30175,31465,31864,33001,35619,37413]
MD=pickle.load(open('full_model.pkl','rb')); NODE=MD['NODE']; OUT=MD['OUT']; ROOT=MD['ROOT']; tree=MD['tree']
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
sys.setrecursionlimit(100000)
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
freeall={}
def fa(v):
    if v in freeall: return freeall[v]
    if v not in defrhs: freeall[v]={v}; return freeall[v]
    freeall[v]=set(); s=set()
    for u in vars_of(defrhs[v]): s|=fa(u)
    freeall[v]=s; return s
atomu={}
for a in E.res:
    s=set()
    for v in vars_of(E.atoms[a]): s|=fa(v)
    hs=[v for v in s if v in handle]
    if len(hs)==1: atomu[a]=hs[0]
u2atom={u:a for a,u in atomu.items()}
h2atom={}
for a,u in atomu.items():
    for v in set(vars_of(E.atoms[a])):
        if v in defrhs and fa(v)=={u}: h2atom.setdefault(v,a)
# selector / gate maps
sel2node={}
for n in NODE:
    for k in ('sa','sb','sab'):
        for s in NODE[n][k]: sel2node[s]=(n,k)
wire2role=collections.defaultdict(list)
for n in NODE:
    for i,d in enumerate(OUT[n]):
        for k in ('va','vb','vab','out'): wire2role[d[k]].append((n,k,i))
# node owning each atom: via any wire in it
def owner(a):
    r=[]
    for v in set(vars_of(E.atoms[a])):
        if v in sel2node: r.append(('gate',)+sel2node[v])
        for t in wire2role.get(v,[]): r.append(('wire',)+t)
    return r
print('%-8s %-9s %s'%('var','kind','detail'))
groups=collections.defaultdict(list)
for v in M32:
    if v in sel2node:
        n,k=sel2node[v]; print('%-8d %-9s selector %s of node x%d'%(v,'GATE',k,n)); groups[n].append((v,'gate:'+k))
    elif v in h2atom:
        a=h2atom[v]; ow=owner(a)
        nodes={t[1] for t in ow}
        print('%-8d %-9s handle of atom %-46s nodes=%s'%(v,'HANDLE',a[:46],sorted(nodes)))
        for nn in nodes: groups[nn].append((v,'handle'))
    elif v in atomu.values():
        a=u2atom[v]; ow=owner(a); nodes={t[1] for t in ow}
        print('%-8d %-9s COFACTOR u of atom %-40s nodes=%s'%(v,'COFACT',a[:40],sorted(nodes)))
        for nn in nodes: groups[nn].append((v,'cofactor'))
    else:
        free = v not in defrhs
        roles=wire2role.get(v,[])
        print('%-8d %-9s free=%s roles=%s def=%s'%(v,'OTHER',free,roles[:2],node_str(defrhs[v])[:50] if not free else '-'))
        for (n,k,i) in roles: groups[n].append((v,'wire:'+k))
print('\nCLUSTERING of the 32 by node:')
for n,vs in sorted(groups.items(),key=lambda kv:-len(kv[1])):
    print('  node x%-6d : %d  %s'%(n,len(set(x[0] for x in vs)),sorted(set(vs))))
un=[v for v in M32 if v not in {x[0] for vs in groups.values() for x in vs}]
print('  unassigned to any node:',un)
