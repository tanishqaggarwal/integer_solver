import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
t=time.time()
sc,v,nz=closure2(BASE,verbose=True)
print('score',sc,'nz',len(nz),'t=%.1f'%(time.time()-t))
for a in nz:
    print('  a%-6d eqs=%-3d %s'%(a,len(L.atom2eq.get(a,{})),L.atom_src[a][:170]))
json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},open('/home/user/integer_solver/solve_lab/agentC_work/c2_%d.json'%sc,'w'))
print('x_11150,x_25739,x_37758 mod p:',v[11150]%P,v[25739]%P,v[37758]%P)
for x in [19083,1308,8778,6418,16742,14853]:
    a=outs.get(x); print('  x_%-6d free=%-5s %s'%(x,x not in outs,(L.atom_src[a][:120] if a is not None else 'FREE')))
