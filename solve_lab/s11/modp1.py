"""Can ANY free variable move x11150, x25739, x37758 modulo p?

If not, the wall congruences are decided entirely by the bit configuration, and no amount of
handle-fiddling in this channel can meet them.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix7_29539_7930.json'))
FREE=[u for u in range(L.NVARS) if u not in L.definer]
T=[11150,25739,37758]
base=[v[u]%P for u in T]
movers={t:[] for t in T}
t0=time.time()
for i,u in enumerate(FREE):
    v2=list(v); L.ripple(v2,{u:v[u]+1})
    for t in T:
        if (v2[t]-v[t])%P: movers[t].append(u)
    if i%1500==0: print(f"  {i}/{len(FREE)} ({time.time()-t0:.0f}s)", flush=True)
for t in T:
    print(f"x{t}: {len(movers[t])} free variables move it mod p -> {movers[t][:20]}")
