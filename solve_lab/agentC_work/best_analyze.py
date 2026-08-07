import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*L.NVARS
for k,val in D.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
av=L.all_atom_values(v)
f=L.failing_eqs(av)
print('score',L.NEQ-len(f),'failing',f)
nzc=[a for a in CHK if av[a]!=0]
nzg=[a for a in L.atom_out if av[a]!=0]
print('nonzero CHECK atoms:',len(nzc),nzc)
print('nonzero GATE atoms:',len(nzg),nzg)
for a in nzc+nzg:
    print(' a%-6d eqs=%-3d out=%-14s val=%s'%(a,len(L.atom2eq.get(a,{})),str(L.atom_out.get(a)),str(av[a])[:40]))
    print('    ',L.atom_src[a][:200])
# which equations fail and what atoms they contain
print()
for e in f:
    m,sq,co=L.eq_atoms[e]
    print('eq %d sq=%s atoms(nonzero only)=%s'%(e,sq,{a:c for a,c in co.items() if av[a]!=0}))
# key selector state
for x in [7715,34554,15298,542,91,22162,30213,16742,14853,12186,24908]:
    print('  x_%-6d = %s'%(x,str(v[x])[:40]))
