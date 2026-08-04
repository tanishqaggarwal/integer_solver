#!/usr/bin/env python3
"""Sweep activators within a quadrant. Records satisfied, x_15298, S0/T0, core/noncore fails.
Usage: python3 agentD_scan_act.py <quadrant 10|01|11> [limit]
Saves any config with satisfied>39016 to best_agentD_<sat>.json and dumps full table json."""
import sys, json, time
import agentD_harness as H
C1,C2=H.C1,H.C2
quad=sys.argv[1] if len(sys.argv)>1 else '10'
limit=int(sys.argv[2]) if len(sys.argv)>2 else 10**9
CONST={30213:C2, 22162:C1, 24468:C1, 18956:C2}
act7=json.load(open('act7715.json'))['free7']
act34=json.load(open('act34554.json'))
rows=[]
best=(39016,None,None)  # (satisfied, override, val)
t0=time.time()
def do(tag, ov):
    global best
    r=H.run_config(ov)
    row={'tag':tag,'sat':r['satisfied'],'nfail':r['nfail'],'x15298':r['x_15298'],
         'core_fail':r['core_fail'],'noncore':r['noncore_fail'],'S0':r['S_is0'],'T0':r['T_is0'],
         'Smodp':str(r['S_modp']),'Tmodp':str(r['T_modp']),'F':r['F']}
    rows.append(row)
    if r['satisfied']>best[0]:
        r2=H.run_config(ov, want_val=True)
        best=(r['satisfied'], ov, r2['val'])
        json.dump({f"x_{i}":r2['val'][i] for i in range(H.NVARS)}, open(f"best_agentD_{r['satisfied']}.json",'w'))
        print(f"  *** NEW BEST {r['satisfied']} at {tag} -> saved", flush=True)
    return r

if quad=='10':
    cands=act7[:limit]
    for i,a in enumerate(cands):
        r=do(f"10:a7={a}", {a:1, **CONST})
        if i%10==0 or r['noncore_fail']<3 or r['satisfied']>39016:
            print(f"[{i}/{len(cands)}] a7={a}: sat={r['satisfied']} core={r['core_fail']} noncore={r['noncore_fail']} S0={r['S_is0']} T0={r['T_is0']} [{time.time()-t0:.0f}s]", flush=True)
elif quad=='01':
    cands=act34[:limit]
    for i,a in enumerate(cands):
        r=do(f"01:a34={a}", {a:1, **CONST})
        if i%10==0 or r['noncore_fail']<3 or r['satisfied']>39016:
            print(f"[{i}/{len(cands)}] a34={a}: sat={r['satisfied']} core={r['core_fail']} noncore={r['noncore_fail']} S0={r['S_is0']} T0={r['T_is0']} [{time.time()-t0:.0f}s]", flush=True)
elif quad=='11':
    # sample pairs
    a7s=act7[:limit]; a34s=act34[:limit]
    for i,a in enumerate(a7s):
        for j,b in enumerate(a34s):
            r=do(f"11:a7={a},a34={b}", {a:1, b:1, **CONST})
            if r['satisfied']>39013:
                print(f"a7={a} a34={b}: sat={r['satisfied']} core={r['core_fail']} noncore={r['noncore_fail']}", flush=True)
        print(f"[{i}/{len(a7s)}] a7={a} done [{time.time()-t0:.0f}s]", flush=True)

# summary
bycount={}
for row in rows:
    bycount.setdefault(row['sat'],0); bycount[row['sat']]+=1
print("\nDISTRIBUTION of satisfied counts:", dict(sorted(bycount.items(), reverse=True)), flush=True)
# best rows
rows.sort(key=lambda r:-r['sat'])
print("TOP 10 configs:")
for row in rows[:10]:
    print(f"  {row['tag']}: sat={row['sat']} core={row['core_fail']} noncore={row['noncore']} S0={row['S0']} T0={row['T0']}", flush=True)
json.dump(rows, open(f"agentD_scan_{quad}.json",'w'))
print(f"\nBEST: {best[0]}  saved table agentD_scan_{quad}.json  [{time.time()-t0:.0f}s]", flush=True)
