import sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
v=forward([0]*L.NVARS)
seen=set()
def show(x,d=0,maxd=8):
    if d>maxd: print('  '*d+'x_%d ...'%x); return
    a=outs.get(x)
    if a is None:
        print('  '*d+'x_%d = FREE (val %s)'%(x,v[x])); return
    src=L.atom_src[a]
    print('  '*d+'x_%d = [a%d] %s   (val %s)'%(x,a,src[:160],v[x]))
    if x in seen: print('  '*d+'  (seen)'); return
    seen.add(x)
    for u in sorted(L.avars[a]):
        if u!=x: show(u,d+1,maxd)
for root in [int(a) for a in sys.argv[1:]]:
    seen=set(); print('======== ROOT x_%d'%root); show(root)
