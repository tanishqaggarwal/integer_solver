"""Full closure scan: one bit at a time, several freeze policies; report best score."""
import ev, fast, json, os, time, collections, sys
from fast import St
from close2 import close
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
st0=St({})
POL={
 'plain':      (set(ALL), {}),
 'pinres':     (set(ALL)|{18956,24468}, {18956:C1,24468:C2}),
}
best={}
res=[]
t0=time.time()
for name,(fr,extra) in POL.items():
    for side in ('A','B'):
        for b in BITS[side]:
            ch=dict(extra); ch[b]=1
            st=st0.clone().set_free(ch)
            out,ok,_=close(st,frozen=set(fr),maxsteps=300)
            res.append((out.score(),name,side,b,tuple(sorted(out.nz()))))
res.sort(reverse=True)
for r in res[:20]: print(r)
print(collections.Counter((r[0],r[1]) for r in res).most_common(10))
print('%.1fs'%(time.time()-t0))
json.dump([[r[0],r[1],r[2],r[3],list(r[4])] for r in res],open('scan1.json','w'))
