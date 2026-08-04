#!/usr/bin/env python3
import json,re
from collections import defaultdict,deque,Counter
p=2**256-2**32-977; NVARS=38748
VAR=re.compile(r'x_(\d+)')
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'],d['rhs'],tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=open('../EQUATIONS.txt').read().split('\n')
broken35=[181,1716,2358,3350,5324,7896,9041,9763,9911,11226,14125,14673,14879,15558,17257,19997,21000,21237,22341,22534,22997,24177,24707,25166,28635,28929,29330,32026,33409,33787,35205,35512,36760,37701,38051]
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
# For the 35, which free inputs (besides 4432,7068) appear?
freecnt=Counter()
for i in broken35:
    vs=set(int(m) for m in VAR.findall(lines[i]))
    for v in vs:
        if v in freeinp and v not in (4432,7068): freecnt[v]+=1
print(f"free inputs (other than 4432,7068) in the 35 broken eqs:")
for v,c in freecnt.most_common(30): print(f"  x_{v}: in {c} of the 35")
print(f"total distinct such free inputs: {len(freecnt)}")
# Also: do the 35 each contain at least one such free knob?
nofree=[i for i in broken35 if not (set(int(m) for m in VAR.findall(lines[i]))&freeinp - {4432,7068})]
print(f"\nof the 35, how many have NO other free input: {len(nofree)} -> {nofree}")
