import sys, os, json, re, collections
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
d=json.load(open('ortree_8599.json'))
for v in d['leaves'][:6]:
    print('LEAF x%d  defined=%s'%(v, v in defrhs))
    for a in resby.get(v,[]):
        print('   ', a)
    print()
