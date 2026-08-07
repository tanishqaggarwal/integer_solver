#!/usr/bin/env python3
import sys,os,json,pickle,time,random
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
import gs2
E=gs2.E
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
K2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
tgt=[i for i,a in enumerate(E.res) if 'x11150' in a or 'x25739' in a or 'x37758' in a]
X3,Y3=22162,30213
def required(base):
    def val(assign):
        v=list(base)
        for k,x in assign.items(): v[k]=x
        r=E.run(v); return [r[i]%p for i in tgt]
    b0=val({})
    cx=[(val({X3:base[X3]+1})[k]-b0[k])%p for k in range(3)]
    cy=[(val({Y3:base[Y3]+1})[k]-b0[k])%p for k in range(3)]
    det=(cx[0]*cy[1]-cx[1]*cy[0])%p
    if det==0: return None
    di=pow(det,p-2,p)
    dx=((-b0[0])*cy[1]+b0[1]*cy[0])%p*di%p
    dy=(cx[0]*(-b0[1])+cx[1]*b0[0])%p*di%p
    if (cx[2]*dx+cy[2]*dy+b0[2])%p: return 'rank3'
    return (base[X3]+dx)%p,(base[Y3]+dy)%p
def run_pair(ba,bb,verbose=False):
    v=[0]*NV
    for k,x in {ba:1,bb:1,22162:K1,30213:K2,24468:K1,18956:K2}.items(): v[k]=x
    v,ok=gs2.solve(v,verbose=verbose,frozen={22162,30213,24468,18956,ba,bb})
    r=E.run(v); bad=E.score(r)
    nz=[i for i in range(len(r)) if r[i]]
    return v,39033-len(bad),nz
if __name__=='__main__':
    pins=json.load(open(os.path.join(HERE,'pins.json')))
    sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
    A=[b for b in sup['7715'] if str(b) in pins]; B=[b for b in sup['34554'] if str(b) in pins]
    random.seed(11)
    out=[]
    for t in range(8):
        ba=random.choice(A); bb=random.choice(B)
        t0=time.time()
        v,s,nz=run_pair(ba,bb)
        rq=required(v)
        rec=dict(ba=ba,bb=bb,score=s,nnz=len(nz),
                 x1=v[12186]%p,y1=v[16742]%p,x2=v[14853]%p,y2=v[24908]%p,
                 pa=[c for _,c in pins[str(ba)]],pb=[c for _,c in pins[str(bb)]],
                 req=rq,t=round(time.time()-t0,1))
        print(json.dumps({k:(str(x) if isinstance(x,int) else x) for k,x in rec.items()})[:600],flush=True)
        out.append(rec)
    pickle.dump(out,open(os.path.join(HERE,'pairs.pkl'),'wb'))
