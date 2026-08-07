import sys, os, json, re, collections, pickle
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
D=pickle.load(open('ortree2.pkl','rb'))
tree=D['tree']
NODE=pickle.load(open('nodes.pkl','rb'))
OUT=pickle.load(open('outwires.pkl','rb'))
def isfree(v): return v not in defrhs
cands=[n for n,N in NODE.items() if tree[N['a']] is None and tree[N['b']] is None and isfree(N['a']) and isfree(N['b'])]
print('bottom nodes with BOTH children free leaves:',len(cands))
n=cands[0]; N=NODE[n]
print('NODE x%d  a=x%d  b=x%d'%(n,N['a'],N['b']))
print(' sel_a',N['sa'],'sel_b',N['sb'],'sel_ab',N['sab'])
for k in ('a','b'):
    print(' leaf x%d atoms:'%N[k])
    for a in resby.get(N[k],[]): print('    ',a[:230])
for d in OUT[n]:
    print(' coord: va=x%d vb=x%d vab=x%d out=x%d'%(d['va'],d['vb'],d['vab'],d['out']))
    for w in (d['va'],d['vb'],d['vab']):
        print('   wire x%d free=%s'%(w,isfree(w)))
        for a in resby.get(w,[]): print('        ',a[:230])
