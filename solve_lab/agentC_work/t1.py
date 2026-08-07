import time,sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
t=time.time(); v=forward([0]*L.NVARS); print('fwd',time.time()-t)
t=time.time(); s,av=score(v); print('score',s,time.time()-t)
# find definer of x_9274 and inspect the 3 core conditions' variables
for x in [9274,29237,23134,7715,34554,18956,14257,24468,32989]:
    a=outs.get(x)
    print('x_%-6d val=%-25s definer a%s : %s'%(x,v[x],a,(L.atom_src[a][:150] if a is not None else 'FREE')))
