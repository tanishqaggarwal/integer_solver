import sys, re, json, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626%P
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002%P
lv=set()
for r in [8599,21839,25956,7304]: lv|=set(leaves(r))
fr=sorted(x for x in lv if x not in outs)
NUM=re.compile(r'x_(\d+) \* \(x_(\d+) - (-?\d+)\)')
info={}
for b in fr:
    pins=[]
    for a in L.var_atoms[b]:
        if a in L.atom_out: continue
        m=NUM.match(L.atom_src[a])
        if m and int(m.group(1))==b: pins.append((int(m.group(2)),int(m.group(3))%P))
    info[b]=pins
two=[b for b in fr if len(info[b])==2]
print('bits with exactly 2 pins:',len(two))
def fit(p1,p2):
    (x1,y1),(x2,y2)=p1,p2
    # y^2 - x^3 = a*x + b
    r1=(y1*y1-x1*x1*x1)%P; r2=(y2*y2-x2*x2*x2)%P
    if (x1-x2)%P==0: return None
    a=(r1-r2)*pow(x1-x2,P-2,P)%P
    b=(r1-a*x1)%P
    return a,b
def on(a,b,x,y): return (y*y-x*x*x-a*x-b)%P==0
cands={}
b0,b1=two[0],two[1]
c0=[c for _,c in info[b0]]; c1=[c for _,c in info[b1]]
for o0 in [(0,1),(1,0)]:
    for o1 in [(0,1),(1,0)]:
        r=fit((c0[o0[0]],c0[o0[1]]),(c1[o1[0]],c1[o1[1]]))
        if not r: continue
        a,b=r
        # test on more bits
        good=0
        for bb in two[2:20]:
            cc=[c for _,c in info[bb]]
            if on(a,b,cc[0],cc[1]) or on(a,b,cc[1],cc[0]): good+=1
        print('order',o0,o1,'a=',a,'b=',b,'matches',good,'/18','  Q on curve?',on(a,b,K2,K1),on(a,b,K1,K2))
