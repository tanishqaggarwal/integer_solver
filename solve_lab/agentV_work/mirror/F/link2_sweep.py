#!/usr/bin/env python3
"""Exhaustive check of link (ii): for EVERY boolean j, turning it on forces the selected wire pair
to j's two pin constants mod p.  Repairs the chains, then compares."""
import sys,os,json,pickle,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import gs2
from fwd import NV
E=gs2.E
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
K2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
pins=json.load(open(os.path.join(HERE,'pins.json')))
sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
A=[b for b in sup['7715'] if str(b) in pins]; B=[b for b in sup['34554'] if str(b) in pins]
outp=os.path.join(HERE,'link2_results.json')
res=json.load(open(outp)) if os.path.exists(outp) else {}
def run(bj, tree):
    partner = 5090 if tree=='A' else 22106
    v=[0]*NV
    for k,x in {bj:1,partner:1,22162:K1,30213:K2,24468:K1,18956:K2}.items(): v[k]=x
    v,ok=gs2.solve(v,verbose=False,frozen={22162,30213,24468,18956,bj,partner})
    if tree=='A': cx,cy=v[12186]%p, v[16742]%p
    else:         cx,cy=v[14853]%p, v[24908]%p
    cs={c%p for _,c in pins[str(bj)]}
    return {'tree':tree,'coord':[str(cx),str(cy)],'match':(cx in cs and cy in cs),'ncs':len(cs)}
todo=[(b,'A') for b in A]+[(b,'B') for b in B]
t0=time.time()
for bj,tree in todo:
    if str(bj) in res: continue
    try: r=run(bj,tree)
    except Exception as ex: r={'tree':tree,'error':str(ex)}
    res[str(bj)]=r
    json.dump(res,open(outp,'w'))
    n=len(res); good=sum(1 for x in res.values() if x.get('match'))
    print('%4d/%d  bit %-6s tree %s match=%s   good=%d  t=%.0fs'%(n,len(todo),bj,tree,r.get('match'),good,time.time()-t0),flush=True)
print('DONE: %d of %d booleans confirm link (ii)'%(sum(1 for x in res.values() if x.get('match')),len(res)))
