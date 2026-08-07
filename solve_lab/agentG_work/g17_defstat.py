import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
v=L.load('/home/user/integer_solver/solve_lab/s10/AG_39013.json'); ad.fwd(v,rounds=6)
SYMS=set(json.load(open('closed_nonbool.json')))
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
def isbool(u):
    for a in L.var_atoms[u]:
        pl=L.polys[a]
        if len(pl)==2 and (u,) in pl and (u,u) in pl and pl[(u,)]==-pl[(u,u)]: return True
    return False
for u in [14257,32989,29309,26777,13458,36358,9899,29967,35795,9254,7927,36864,25295,18956,24468,2081,24601,12378,30033]:
    st = 'DEFINED by a%d'%L.definer[u] if u in L.definer else ('FREE'+(' BOOL' if isbool(u) else '')+(' [in SYMS]' if u in SYMS else ' [NOT in SYMS]'))
    print('x%-6d natoms=%-3d val=%s  %s'%(u,len(L.var_atoms[u]),str(v[u])[:16],st))
