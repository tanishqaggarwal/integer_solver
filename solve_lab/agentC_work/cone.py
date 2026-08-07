import sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
v=forward([0]*L.NVARS)
def cone(roots,maxn=100000):
    seen=set(); st=list(roots)
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        a=outs.get(x)
        if a is None: continue
        for u in L.avars[a]:
            if u!=x and u not in seen: st.append(u)
    return seen
roots=[int(a) for a in sys.argv[1:]]
c=cone(roots)
print('cone size',len(c),'free in cone',len([x for x in c if x not in outs]))
tp={u:i for i,u in enumerate(topo)}
for x in sorted(c,key=lambda z:tp[z]):
    a=outs.get(x)
    if a is None: print('x_%-6d = FREE   val=%s'%(x,v[x]))
    else: print('x_%-6d = %s   | val=%s'%(x,L.atom_src[a][:170],v[x]))
