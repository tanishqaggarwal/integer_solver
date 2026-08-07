import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
from ec import add,neg,mul,leafpoints,N
from ort import leaves
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
chain=[int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentC_work/chain.json'))]
pts=leafpoints()
side={}
for r,tag in [(8599,'s1'),(21839,'s1'),(25956,'s2'),(7304,'s2')]:
    for x in leaves(r):
        if x in pts: side[x]=tag
bitsA=[e for e in range(256) if (N>>e)&1 and side[chain[e]]=='s1']
bitsB=[e for e in range(256) if (N>>e)&1 and side[chain[e]]=='s2']
print('popcount(n)=',bin(N).count('1'),'on s1:',len(bitsA),'on s2:',len(bitsB))
kA=sum(1<<e for e in bitsA); kB=sum(1<<e for e in bitsB)
PA=mul(kA,pts[chain[0]]); PB=mul(kB,pts[chain[0]])
print('PA+PB == O ?', add(PA,PB) is None)
print('PA',PA)
ctrl={22162:K2, 30213:K1}
for e in bitsA+bitsB: ctrl[chain[e]]=1
t=time.time()
sc,v,nz=closure2(ctrl,rounds=25,verbose=True)
print('phase1 score',sc,'nz',len(nz),'t=%.1f'%(time.time()-t))
print('x1=%d'%(v[12186]%P)); print('x2(x_1308)=%d'%(v[1308]%P))
print('y1(x_19083)=%d'%(v[19083]%P)); print('y2(x_24908)=%d'%(v[24908]%P))
print('PA.x==x_12186?',PA[0]==v[12186]%P,' PA.y==x_19083?',PA[1]==v[19083]%P)
for a in nz[:20]:
    print('  a%-6d eqs=%-3d %s'%(a,len(L.atom2eq.get(a,{})),L.atom_src[a][:150]))
json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},open('/home/user/integer_solver/solve_lab/agentC_work/np_phase1_%d.json'%sc,'w'))
