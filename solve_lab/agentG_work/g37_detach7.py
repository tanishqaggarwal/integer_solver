import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gframe
from gsym2 import L, ad, P
A7=[22229,22230,35758,35759,35760,35761,35762]
for a in A7:
    print('a%-6d gate=%s out=%s neq=%d'%(a,a in L.atom_out,L.atom_out.get(a),len(L.atom2eq.get(a,{}))))
D=[L.atom_out[a][1] for a in A7 if a in L.atom_out]
print('detach targets',D)
v=L.load('/home/user/integer_solver/solve_lab/s10/AG_39013.json'); ad.fwd(v,rounds=6)
nfail,inc,res=gframe.ceiling(v,D)
print('ceiling with those detached: %d (fail %d)'%(L.NEQ-nfail,nfail))
print('residual equations:',[i for i,_ in res][:30],'inc',inc[:10])
