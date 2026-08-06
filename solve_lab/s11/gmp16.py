"""Base with the two free checks already cleared, then characterise the 256 REAL bits.

1,156 free bits: 900 are inert (flipping them changes no residue at all), so the genuine search
space is the other 256 -- the message bits.  x9118 kills a35759 and x8731 kills a35760 with no
collateral whatsoever, so start from the four-failing state and score every real bit there.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
BITS=[]
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.append(ks[0][0])
BITS=sorted(set(u for u in BITS if u in FREE))
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
forwardp(base)
def fails(v): return [a for a in CHK if evalp(L.polys[a],v)]
# clear a35759 via x9118 and a35760 via x8731 (measured responses: 5113045 and -1)
for a,u,c in [(35759,9118,5113045),(35760,8731,P-1)]:
    r=evalp(L.polys[a],base)
    base[u]=(base[u]-r*pow(c,-1,P))%P
    forwardp(base)
F=fails(base)
print("after clearing the two free checks:",len(F),F)
json.dump([int(x) for x in base], open(os.path.join(HERE,'data','gmp16_base.json'),'w'))
inert=set(u for u,*_ in [tuple(x) for x in json.load(open(os.path.join(HERE,'data','gmp15.json')))])
real=[u for u in BITS if u not in inert]
print(f"real (non-inert) bits: {len(real)}   currently ON: {[u for u in real if base[u]%P==1]}")
t0=time.time(); res=[]
for u in real:
    v=list(base); v[u]=(1-base[u])%P; forwardp(v)
    f=fails(v)
    res.append((len(f),u,base[u]%P,f[:6]))
res.sort()
print(f"scanned {len(res)} real bits ({time.time()-t0:.0f}s)")
print("  distribution:", dict(sorted(collections.Counter(r[0] for r in res).items())))
for n,u,old,f in res[:12]:
    print(f"   x{u}: {old}->{1-old}  failing={n}  {f}")
json.dump([[r[0],r[1],r[2]] for r in res], open(os.path.join(HERE,'data','gmp16.json'),'w'))
