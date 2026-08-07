import sys, json, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close import closure, evalat, CHK
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
def state(extra):
    s=dict(BASE); s.update(extra)
    seeds,v,sc,nz=closure(s,verbose=False)
    return seeds,v,sc,nz
def rep(v,tag):
    x1,y1,x2,y2,x3,y3=[v[i]%P for i in (12186,16742,14853,24908,22162,30213)]
    A=v[35389]%P; B=v[6671]%P
    Ap=((x2-x1)**2*(x3+x1+x2)-(y2-y1)**2)%P
    Bp=((y3+y1)*(x2-x1)-(x1-x3)*(y2-y1))%P
    print('[%s] A=%d Ap=%d match=%s | B=%d Bp=%d match=%s'%(tag,A,Ap,A==Ap,B,Bp,B==Bp))
    return (x1,y1,x2,y2,x3,y3,A,B)
seeds,v,sc,nz=state({})
print('score',sc,'nz',len(nz))
c=rep(v,'default')
print('coords',[hex(z)[:18] for z in c])
random.seed(7)
for t in range(3):
    r1=random.randrange(P); r2=random.randrange(P)
    seeds,v,sc,nz=state({16742:r1,14853:r2})
    print('  probe',t,'score',sc,'nz',len(nz),'y1 kept',v[16742]%P==r1,'x2 kept',v[14853]%P==r2)
    rep(v,'p%d'%t)
