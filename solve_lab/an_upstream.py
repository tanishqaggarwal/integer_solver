#!/usr/bin/env python3
import heal_harness as H
import json
p=H.p
gates={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.setdefault(d['t'],[]).append((d['rhs'],tuple(d['vids'])))
# fan-out (eqs touched) per var
from collections import defaultdict
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for wv in H.eqvars[i]:
        if wv in H.freeinp: s.add(wv)
        s|=H.anc.get(wv,set())
    eq_free.append(s)
# descendants-of-var -> eqs
def eqs_touched_by(t):
    desc=set([t])
    for tt in H.order:
        if t in H.anc.get(tt,set()): desc.add(tt)
    return set(i for i in range(len(H.eqcode)) if H.eqvars[i]&desc), desc
F16=set([697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431])
F11=set([2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125])
# trace: for x_19964 (=x_20492+x_19892) and x_2099 (=x_37158+x_25297), examine each input's touched eqs
for t,inps in [(19964,[20492,19892]),(2099,[37158,25297])]:
    print(f"=== x_{t} inputs ===")
    for w in inps:
        te,desc=eqs_touched_by(w)
        print(f"  x_{w}: free={w in H.freeinp}, touches {len(te)} eqs; of16={len(te&F16)}, of11={len(te&F11)}, other={len(te-F16-F11)}")
        if w in gates:
            print(f"     def: {gates[w][0][0][:80]}")
