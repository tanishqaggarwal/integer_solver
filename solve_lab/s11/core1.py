"""Which free variables move the core trio without disturbing the mod-p work?"""
import sys, os, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','modp5_out.json'))
FREE=[u for u in range(L.NVARS) if u not in L.definer]
MODP=[11150,25739,37758]
TRIO=[688,1618,40608]
def av(v,a): return L.evalpoly(L.polys[a],v)
base={a:av(v,a) for a in TRIO}
res=[]
t0=time.time()
for u in FREE:
    v2=list(v); L.ripple(v2,{u:v[u]+1})
    d=[av(v2,a)-base[a] for a in TRIO]
    if not any(d): continue
    keeps = all(v2[t]%P==0 for t in MODP)
    res.append((u,keeps,[x!=0 for x in d]))
print(f"free vars moving the core trio: {len(res)}  ({time.time()-t0:.0f}s)")
keep=[r for r in res if r[1]]
print(f"  of which PRESERVE x11150,x25739,x37758 = 0 mod p: {len(keep)}")
for r in keep[:25]: print(f"     x{r[0]} touches {['a688','a1618','a40608'][0] if 0 else ''}", r[2])
