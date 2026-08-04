#!/usr/bin/env python3
"""Sweep activators combined with the WINNING routing (30213=0,22162=0,24468=C1,18956=C2).
Both quadrants (1,0) via act7715 and (0,1) via act34554. Looks for >39018."""
import json, time
import agentD_harness as H
C1,C2=H.C1,H.C2
WIN={30213:0,22162:0,24468:C1,18956:C2}
act7=json.load(open('act7715.json'))['free7']
act34=json.load(open('act34554.json'))
best=(39018,None)
t0=time.time()
dist={}
for tag,acts in [('10',act7),('01',act34)]:
    for i,a in enumerate(acts):
        r=H.run_config({a:1, **WIN})
        dist[r['satisfied']]=dist.get(r['satisfied'],0)+1
        if r['satisfied']>best[0]:
            best=(r['satisfied'],(tag,a))
            r2=H.run_config({a:1,**WIN}, want_val=True)
            json.dump({f"x_{i2}":r2['val'][i2] for i2 in range(H.NVARS)}, open(f"best_agentD_{r['satisfied']}.json",'w'))
            print(f"  *** NEW BEST {r['satisfied']} at {tag} a={a}", flush=True)
        elif r['satisfied']==39018 and i%25==0:
            print(f"  [{tag} {i}] a={a}: 39018 (S0={r['S_is0']} T0={r['T_is0']}) [{time.time()-t0:.0f}s]", flush=True)
print(f"\nDISTRIBUTION: {dict(sorted(dist.items(),reverse=True))}")
print(f"BEST with winning routing: {best}  [{time.time()-t0:.0f}s]")
