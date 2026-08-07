import sys, os, json, re, collections, pickle
F='/home/user/integer_solver/solve_lab/agentV_work/mirror/F'; sys.path.insert(0,F)
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
# pick a node whose child is an OR node
for n,N in NODE.items():
    if tree[N['a']] is not None:
        print('node x%d child_a x%d'%(n,N['a']))
        print(' n reads', [d['va'] for d in OUT[n]])
        print(' child outs', [d['out'] for d in OUT[N['a']]])
        for w in [d['va'] for d in OUT[n]]:
            print('  x%d def=%s'%(w, node_str(defrhs[w])[:80] if w in defrhs else 'FREE'))
            for a in resby.get(w,[]): print('     ',a[:200])
        for o in [d['out'] for d in OUT[N['a']]]:
            print('  child out x%d used in %s'%(o,[(u,node_str(defrhs[u])[:60]) for u in [w for w,r in defrhs.items() if o in vars_of(r)]][:4]))
            for a in resby.get(o,[]): print('     res:',a[:200])
        break
