import sys, itertools, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
CHK=[a for a in range(L.NA) if a not in L.atom_out]
def run(seeds,tag):
    v=[0]*L.NVARS
    for k,val in seeds.items(): v[k]=val
    forward(v)
    for k,val in seeds.items():
        if v[k]!=val: print('   !! seed x_%d overwritten'%k)
    av=L.all_atom_values(v)
    f=L.failing_eqs(av)
    nz=[a for a in CHK if av[a]!=0]
    print('[%s] score=%d failing=%d nzchecks=%d  s1=%s s2=%s  R1=%s R2=%s'%(
        tag,L.NEQ-len(f),len(f),len(nz),v[7715],v[34554],(v[18956]-K1)%P==0,(v[24468]-K2)%P==0))
    return v,av,nz,f
# baseline
run({}, 'base')
run({542:1}, 'x542=1')
run({91:1}, 'x91=1')
run({542:1,91:1}, 'both bits')
v,av,nz,f=run({542:1,91:1,22162:K2,30213:K1}, 'both+coords')
for a in nz[:40]:
    print('   a%-6d eqs=%-3d val=%s'%(a,len(L.atom2eq.get(a,{})), str(av[a])[:60]))
    print('      ',L.atom_src[a][:200])
