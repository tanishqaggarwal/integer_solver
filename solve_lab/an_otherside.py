#!/usr/bin/env python3
import heal_harness as H
p=H.p
# defining gate of x_19964 and x_2099
gd={t:(rhs,vids) for t,rhs,vids in [(g[0],g[1],g[2]) for g in [(d['t'],d['rhs'],d['vids']) for d in __import__('json').loads('[]')]]} if False else None
import json
gates={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.setdefault(d['t'],[]).append((d['rhs'],tuple(d['vids'])))
for t in [19964,2099]:
    print(f"=== x_{t} ===")
    print(f"  free ancestors: {len(H.anc.get(t,set()))} -> {sorted(H.anc.get(t,set()))[:12]}")
    if t in gates:
        for rhs,vids in gates[t][:2]:
            print(f"  gate rhs: {rhs[:150]}")
# Does x_19964 / x_2099 feed the 16 ripple eqs? (are they ancestors of the 16's vars)
F16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
# compute descendants of x_19964, x_2099 among gate outputs, then eqs
for t in [19964,2099]:
    desc=set([t])
    for tt in H.order:
        if t in H.anc.get(tt,set()): desc.add(tt)
    hit=[i for i in F16 if H.eqvars[i]&desc]
    print(f"\nx_{t} feeds {len(desc)} vars; touches {len(hit)} of the 16 ripple eqs: {hit}")
    # also touches the 11?
    F11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
    hit11=[i for i in F11 if H.eqvars[i]&desc]
    print(f"    touches {len(hit11)} of the 11 obstruction eqs")
