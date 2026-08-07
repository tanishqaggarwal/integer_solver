import sys, json, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
sc,v,nz=closure2(BASE)
KA=97553848499418123410591666447050222001188385549510401465815187079080512838891
x1=v[12186]%P; y1=v[16742]%P; x2=v[14853]%P; y2=v[24908]%P; x3=v[22162]%P; y3=v[30213]%P
pts=[(x1,y1),(x2,y2),(x3,y3)]
print('P1',x1,y1); print('P2',x2,y2); print('P3',x3,y3)
# y^2 = x^3 + KA x^2 + a4 x + a6
def rhs0(x): return (pow(x,3,P)+KA*x*x)%P
r1=(y1*y1-rhs0(x1))%P; r2=(y2*y2-rhs0(x2))%P; r3=(y3*y3-rhs0(x3))%P
a4=(r1-r2)*pow(x1-x2,P-2,P)%P; a6=(r1-a4*x1)%P
print('a4',a4); print('a6',a6)
print('P3 on curve?',(r3-a4*x3-a6)%P==0)
# check P1+P2 = P3 with this law
lam=(y2-y1)*pow(x2-x1,P-2,P)%P
X3=(lam*lam-KA-x1-x2)%P; Y3=(lam*(x1-X3)-y1)%P
print('addition gives x3=%s (match %s), y3 match %s'%(X3,X3==x3,Y3==y3))
json.dump({'KA':str(KA),'a4':str(a4),'a6':str(a6),'P1':[str(x1),str(y1)],'P2':[str(x2),str(y2)],'Q':[str(x3),str(y3)]},
          open('/home/user/integer_solver/solve_lab/agentC_work/curve.json','w'))
