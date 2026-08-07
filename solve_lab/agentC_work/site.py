import sys, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
sc,v,nz=closure2(BASE)
print('base',sc)
def show(x,d=0,md=4,seen=None):
    if seen is None: seen=set()
    if x in seen or d>md: return
    seen.add(x)
    a=outs.get(x)
    if a is None:
        print('  '*d+'x_%-6d FREE   val=%s'%(x,str(v[x])[:24])); return
    print('  '*d+'x_%-6d <- a%-6d eqs=%-3d  %s  val=%s'%(x,a,len(L.atom2eq.get(a,{})),L.atom_src[a][:110],str(v[x])[:22]))
    for u in sorted(L.avars[a]):
        if u!=x: show(u,d+1,md,seen)
for a in [26731,29539]:
    print('==== a%d eqs=%d : %s'%(a,len(L.atom2eq.get(a,{})),L.atom_src[a]))
    for u in sorted(L.avars[a]): show(u,1,3)
