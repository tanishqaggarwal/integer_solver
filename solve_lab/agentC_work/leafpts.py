import sys, re, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
print('Q=(K2 mod p, K1 mod p) on y^2=x^3+7 ?', pow(K1%P,2,P)==(pow(K2%P,3,P)+7)%P)
lv=set()
trees={}
for r in [8599,21839,25956,7304]:
    trees[r]=leaves(r); lv|=set(trees[r])
fr=sorted(x for x in lv if x not in outs)
NUM=re.compile(r'x_(\d+) \* \(x_(\d+) - (-?\d+)\)')
info={}
for b in fr:
    pins=[]
    for a in L.var_atoms[b]:
        if a in L.atom_out: continue
        s=L.atom_src[a]
        m=NUM.match(s)
        if m and int(m.group(1))==b:
            pins.append((a,int(m.group(2)),int(m.group(3)),s))
    info[b]=pins
cnt=collections.Counter(len(v) for v in info.values())
print('pin count hist',cnt)
oncurve=0; pts={}
for b,pins in info.items():
    cs=[c%P for _,_,c,_ in pins]
    found=None
    for i in range(len(cs)):
        for j in range(len(cs)):
            if i==j: continue
            x,y=cs[i],cs[j]
            if pow(y,2,P)==(pow(x,3,P)+7)%P: found=(x,y,pins[i][1],pins[j][1])
    if found: oncurve+=1; pts[b]=found
print('leaf bits whose two pin constants form a secp256k1 point:',oncurve,'of',len(fr))
json.dump({str(k):[str(x) for x in v] for k,v in pts.items()}, open('/home/user/integer_solver/solve_lab/agentC_work/leafpts.json','w'))
# print a few
for b in sorted(pts)[:5]:
    x,y,vx,vy=pts[b]; print(b,'x=',x,'y=',y)
