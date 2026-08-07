import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close3 import *
from pairsweep import build, lp, S1, S2
W='/home/user/integer_solver/solve_lab/agentC_work/'
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
t=time.time()
sc,v,nz=closure3({542:1,91:1,22162:K2,30213:K1},verbose=True)
print('baseline closure3:',sc,'nz',len(nz),'%.1fs'%(time.time()-t))
for a in nz: print('   a%-6d eqs=%-3d %s'%(a,len(L.atom2eq.get(a,{})),L.atom_src[a][:130]))
t=time.time()
sc,v,nz=closure3(build(24601,2081,0),verbose=True)
print('pair 24601/2081:',sc,'nz',len(nz),'%.1fs'%(time.time()-t))
for a in nz: print('   a%-6d eqs=%-3d %s'%(a,len(L.atom2eq.get(a,{})),L.atom_src[a][:130]))
json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},open(W+'c3_%d.json'%sc,'w'))
