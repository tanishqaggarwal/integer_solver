"""Measure the delivered tree point for each single bit, by running the bottom-up closer."""
import ev, fast, json, os, sys, time
from fast import St
from close2 import close
HERE=os.path.dirname(os.path.abspath(__file__))
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
BITS=json.load(open(os.path.join(HERE,'bits.json')))
st0=St({})
res={}
t0=time.time()
for side in ('A','B'):
    xr,yr=(23927,19083) if side=='A' else (1308,17601)
    for b in BITS[side]:
        st=st0.clone().set_free({b:1})
        # freeze the residue-carrying free inputs so the closer cannot fudge them
        out,ok,fr=close(st,frozen={18956,24468,7497,11436,22820,14393,12186,16742,14853,24908,22162,30213},
                        verbose=False,maxsteps=200)
        nz=sorted(out.nz())
        res[b]=(side,out.v[xr]%p,out.v[yr]%p,out.score(),nz)
print('measured %d bits in %.1fs'%(len(res),time.time()-t0))
import collections
print(collections.Counter(tuple(v[4]) for v in res.values()).most_common(6))
json.dump({str(k):[v[0],str(v[1]),str(v[2]),v[3]] for k,v in res.items()},open(os.path.join(HERE,'bitpoints.json'),'w'))
