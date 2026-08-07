import sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
v=forward([0]*L.NVARS)
for x in [int(a) for a in sys.argv[1:]]:
    print('=== x_%d  val=%s  free=%s'%(x,v[x],x not in outs))
    for a in L.var_atoms[x]:
        role='DEF-of-x_%d'%L.atom_out[a][1] if a in L.atom_out else 'CHECK'
        print('   a%-6d %-14s %s'%(a,role,L.atom_src[a][:190]))
