#!/usr/bin/env python3
"""Scan constant-routing combinations for a quadrant. The 4 knobs {30213,22162,24468,18956}
each take a value in {0, C1, C2}. Report satisfied count; save best."""
import sys, json, itertools, time
import agentD_harness as H
C1,C2=H.C1,H.C2
quad=sys.argv[1] if len(sys.argv)>1 else '10'
a7=24601; a34=2081
if quad=='10': base={a7:1}
elif quad=='01': base={a34:1}
elif quad=='11': base={a7:1,a34:1}
elif quad=='00': base={}
knobs=[30213,22162,24468,18956]
vals={'0':0,'C1':C1,'C2':C2}
rows=[]
best=(0,None)
t0=time.time()
combos=list(itertools.product(vals.keys(), repeat=4))
print(f"quad {quad}: scanning {len(combos)} routings", flush=True)
for i,combo in enumerate(combos):
    ov=dict(base)
    for k,c in zip(knobs, combo):
        if vals[c]!=0: ov[k]=vals[c]
        else: ov[k]=0  # explicit 0 pin
    r=H.run_config(ov)
    rows.append((r['satisfied'], combo, r['core_fail'], r['noncore_fail'], r['S_is0'], r['T_is0']))
    if r['satisfied']>best[0]:
        best=(r['satisfied'], combo)
        if r['satisfied']>=39016:
            r2=H.run_config(ov, want_val=True)
            json.dump({f"x_{i2}":r2['val'][i2] for i2 in range(H.NVARS)}, open(f"agentD_route_{quad}_{'_'.join(combo)}.json",'w'))
        print(f"  [{i}] {combo}: sat={r['satisfied']} core={r['core_fail']} noncore={r['noncore_fail']} S0={r['S_is0']} T0={r['T_is0']} *BEST* [{time.time()-t0:.0f}s]", flush=True)
    elif r['satisfied']>=39013:
        print(f"  [{i}] {combo}: sat={r['satisfied']} core={r['core_fail']} noncore={r['noncore_fail']}", flush=True)
rows.sort(reverse=True)
print(f"\nTOP 12 routings for quad {quad}:")
for sat,combo,cf,ncf,s0,t0_ in rows[:12]:
    print(f"  {combo}: sat={sat} core={cf} noncore={ncf} S0={s0} T0={t0_}")
print(f"BEST: {best}  [{time.time()-t0:.0f}s]", flush=True)
json.dump([[s,list(c),cf,ncf,bool(s0),bool(t0_)] for s,c,cf,ncf,s0,t0_ in rows], open(f"agentD_route_{quad}.json",'w'))
