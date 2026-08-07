import sys, re, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
a2=int(C['KA']); a4=int(C['a4']); a6=int(C['a6'])
def on(x,y): return (y*y-(pow(x,3,P)+a2*x*x+a4*x+a6))%P==0
b2=4*a2%P; b4=2*a4%P; b6=4*a6%P; b8=(4*a2*a6-a4*a4)%P
D=(-b2*b2%P*b8 - 8*pow(b4,3,P) - 27*b6*b6 + 9*b2*b4%P*b6)%P
print('discriminant',D,'SINGULAR' if D==0 else 'nonsingular')
c4=(b2*b2-24*b4)%P
print('j =',(pow(c4,3,P)*pow(D,P-2,P))%P if D else None)
def add(Pt,Qt):
    if Pt is None: return Qt
    if Qt is None: return Pt
    x1,y1=Pt; x2,y2=Qt
    if (x1-x2)%P==0:
        if (y1+y2)%P==0: return None
        lam=(3*x1*x1+2*a2*x1+a4)*pow(2*y1,P-2,P)%P
    else:
        lam=(y2-y1)*pow(x2-x1,P-2,P)%P
    x3=(lam*lam-a2-x1-x2)%P; y3=(lam*(x1-x3)-y1)%P
    return (x3,y3)
def mul(k,Pt):
    R=None; Q=Pt
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
P1=(int(C['P1'][0]),int(C['P1'][1])); P2=(int(C['P2'][0]),int(C['P2'][1])); Q=(int(C['Q'][0]),int(C['Q'][1]))
print('P1 on',on(*P1),'P2 on',on(*P2),'Q on',on(*Q))
print('P1+P2 ==Q?',add(P1,P2)==Q)
# leaf points
lv=set()
for r in [8599,21839,25956,7304]: lv|=set(leaves(r))
fr=sorted(x for x in lv if x not in outs)
NUM=re.compile(r'x_(\d+) \* \(x_(\d+) - (-?\d+)\)')
pts={}
for b in fr:
    cs=[]
    for a in L.var_atoms[b]:
        if a in L.atom_out: continue
        m=NUM.match(L.atom_src[a])
        if m and int(m.group(1))==b: cs.append((int(m.group(2)),int(m.group(3))%P))
    found=None
    for i in range(len(cs)):
        for j in range(len(cs)):
            if i!=j and on(cs[i][1],cs[j][1]): found=(cs[i][1],cs[j][1],cs[i][0],cs[j][0])
    if found: pts[b]=found
print('leaf bits on curve:',len(pts),'of',len(fr))
json.dump({str(k):[str(v[0]),str(v[1]),v[2],v[3]] for k,v in pts.items()},open('/home/user/integer_solver/solve_lab/agentC_work/leafpts2.json','w'))
if pts:
    ks=sorted(pts)
    # is it a doubling chain?
    l=[(pts[b][0],pts[b][1]) for b in ks]
    print('first leaf',l[0])
    hits=0
    for i in range(len(l)-1):
        if add(l[i],l[i])==l[i+1]: hits+=1
    print('consecutive-doubling hits (sorted by var id):',hits,'/',len(l)-1)
    # is P1 among leaves?
    print('P1 is a leaf?',P1 in l, 'P2 is a leaf?',P2 in l)
    # anomalous test: is [p]G = O ?
    print('[p]P1 == O ?', mul(P,l[0]) is None)
