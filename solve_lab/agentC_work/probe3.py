import sys, json, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close import closure, evalat, CHK
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
def st(extra):
    s=dict(BASE); s.update(extra)
    seeds,v,sc,nz=closure(s,verbose=False)
    return v,sc,nz
random.seed(11)
rows=[]
for t in range(5):
    e={} if t==0 else {16742:random.randrange(P),14853:random.randrange(P)}
    v,sc,nz=st(e)
    x1,y1,x2,y2,x3,y3=[v[i]%P for i in (12186,16742,14853,24908,22162,30213)]
    A=v[35389]%P; B=v[6671]%P
    d=(x2-x1)%P
    K=((A+(y2-y1)**2)*pow(d*d%P,P-2,P)-(x3+x1+x2))%P if d else None
    print('t=%d sc=%d nz=%d'%(t,sc,len(nz)))
    print('   x1=%d'%x1); print('   y2=%d'%y2); print('   K =%s'%K)
    rows.append((x1,y1,x2,y2,x3,y3,A,B,K))
print()
print('x1 constant across probes:',len(set(r[0] for r in rows))==1)
print('y2 constant across probes:',len(set(r[3] for r in rows))==1)
print('K  constant across probes:',len(set(r[8] for r in rows))==1, set(r[8] for r in rows) if len(set(r[8] for r in rows))<4 else '')
