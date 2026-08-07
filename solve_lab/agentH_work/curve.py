import ev, pickle, json, os
from fast import sup, fidx, inv
HERE=os.path.dirname(os.path.abspath(__file__))
LP=json.load(open(os.path.join(HERE,'leafpins.json')))
BITS=json.load(open(os.path.join(HERE,'bits.json')))
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
# x-side tree roots: A-side x_23927 (->x_12186), B-side x_1308 (->x_14853)
# y-side tree roots: A-side x_19083 (->x_16742), B-side x_17601 (->x_24908)
def has(v,u): return (sup[v]>>fidx[u])&1
pts={}
prob=[]
for b,pins in LP.items():
    b=int(b)
    side='A' if b in BITS['A'] else 'B'
    xr,yr=(23927,19083) if side=='A' else (1308,17601)
    X=Y=None
    for F,H,c in pins:
        if has(xr,F): X=(F,H,c)
        elif has(yr,F): Y=(F,H,c)
    if X is None or Y is None: prob.append(b); continue
    pts[b]=(side,X,Y)
print('classified',len(pts),'problem',prob[:10],len(prob))
json.dump({str(k):v for k,v in pts.items()},open(os.path.join(HERE,'points.json'),'w'))
# curve test in shifted coords x = X + K/3 mod p
inv3=pow(3,p-2,p); sh=K*inv3%p
P=[( (pts[b][1][1]+sh)%p, pts[b][2][1]%p ) for b in sorted(pts)]
# fit y^2 = x^3 + a x + b using first two points
(x1,y1),(x2,y2)=P[0],P[1]
# y1^2 - x1^3 = a x1 + b ; y2^2 - x2^3 = a x2 + b
r1=(y1*y1-pow(x1,3,p))%p; r2=(y2*y2-pow(x2,3,p))%p
a=(r1-r2)*pow((x1-x2)%p,p-2,p)%p
bb=(r1-a*x1)%p
print('fitted a =',a); print('fitted b =',bb)
ok=sum(1 for (x,y) in P if (y*y-pow(x,3,p)-a*x-bb)%p==0)
print('points on that curve:',ok,'of',len(P))
xt=(C2+sh)%p; yt=C1%p
print('target on curve:',(yt*yt-pow(xt,3,p)-a*xt-bb)%p==0)
