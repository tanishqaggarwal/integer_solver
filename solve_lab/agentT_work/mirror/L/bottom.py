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
uses=collections.defaultdict(list)
for w,r in defrhs.items():
    for u in vars_of(r): uses[u].append(w)
D=pickle.load(open('ortree2.pkl','rb'))
tree=D['tree']; selmap=D['selmap']
bottoms=[n for n,ch in tree.items() if ch and tree[ch[0]] is None and tree[ch[1]] is None]
print('bottom OR nodes (both children leaves):',len(bottoms))
n=bottoms[0]; sm=selmap[n]
print('node x%d  a=x%d b=x%d'%(n,sm['a'],sm['b']))
print('  s_a',sm['s_a'],' s_b',sm['s_b'],' s_ab',sm['s_ab'])
print('LEAF a atoms:'); [print('   ',a) for a in resby.get(sm['a'],[])]
print('LEAF b atoms:'); [print('   ',a) for a in resby.get(sm['b'],[])]
for s in sm['s_a']+sm['s_b']+sm['s_ab']:
    print('SEL x%d used in defs of: %s'%(s,uses.get(s)))
    for w in uses.get(s,[]):
        print('    x%d := %s'%(w,node_str(defrhs[w])[:200]))

print()
print('=== x4858 def:', node_str(defrhs.get(4858)) if 4858 in defrhs else 'FREE')
for v in [4843,30800,12398,27914,35518,30176,22551,18213,13496,20593]:
    print('x%d def=%s'%(v, node_str(defrhs[v])[:120] if v in defrhs else 'FREE'))
    print('    used in:', [(w,node_str(defrhs[w])[:100]) for w in uses.get(v,[])][:4])
    print('    res atoms:', [a[:120] for a in resby.get(v,[])][:4])
