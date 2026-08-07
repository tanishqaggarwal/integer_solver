"""Bit-flip scan in the mod-p layer.

A global mod-p forward evaluation costs 0.08s, so every message bit can be flipped and scored
exactly.  This is the first time the bit choice has been measured against the quantity that
actually matters -- the number of checks failing mod p -- rather than against integer score.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
# boolean-checked variables:  atom == x*x - x
BITS=[]
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]):
        u=ks[0][0]
        if Pp[ks[0]]==-Pp[ks[1]]: BITS.append(u)
BITS=sorted(set(BITS))
print("boolean-checked variables:",len(BITS)," free among them:",sum(1 for u in BITS if u in FREE))
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
forwardp(base)
def nf(v): return [a for a in CHK if evalp(L.polys[a],v)]
b0=nf(base)
print("base failing checks mod p:",len(b0),b0)
on=[u for u in BITS if base[u]%P==1]
print("bits currently ON:",len(on),on[:20])
t0=time.time(); res=[]
for u in BITS:
    if u not in FREE: continue
    v=list(base); v[u]=1-base[u]%P if base[u]%P in (0,1) else 1
    v[u]%=P
    forwardp(v)
    f=nf(v)
    res.append((len(f),u,base[u]%P,f[:8]))
res.sort()
print(f"scanned {len(res)} bits ({time.time()-t0:.0f}s).  best flips:")
for n,u,old,f in res[:20]:
    print(f"   x{u}: {old}->{1-old}  failing checks mod p = {n}   {f}")
print("  distribution:", dict(sorted(collections.Counter(r[0] for r in res).items())[:12]))
json.dump([[r[0],r[1],r[2]] for r in res], open(os.path.join(HERE,'data','gmp13.json'),'w'))
