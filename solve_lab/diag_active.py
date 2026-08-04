#!/usr/bin/env python3
"""List the ~30 nonzero free inputs in the best solution. Their values, bit-sizes, and which
control/load role each plays. This is the sparse witness structure."""
import json, re, sys
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
nz=[(v,best[v]) for v in sorted(freeinp) if best.get(v,0)!=0]
print(f"nonzero free inputs: {len(nz)}")
ctrl={14853,12186,16742,30317,5146,2936,24601,2081,30213,22162,24468,18956}
for v,x in nz:
    tag=""
    if v in ctrl: tag=" [CTRL]"
    if x==C1: tag+=" =C1"
    elif x==C2: tag+=" =C2"
    elif x in (1,-1): tag+=" =activator"
    print(f"  x_{v}: {x.bit_length()} bits{tag}")
