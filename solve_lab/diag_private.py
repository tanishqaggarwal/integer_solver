#!/usr/bin/env python3
"""Find handles in the 23 broken eqs that appear in FEW total equations (low ripple). Check if
low-coupling handles span the rank-9 residual space -> clean heal possible."""
import json, re, sys
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
byh=defaultdict(set)
for i,vs in enumerate(eqvars):
    for v in vs&freeinp: byh[v].add(i)
broken=set([3408, 3841, 4134, 4526, 5069, 7276, 15440, 15724, 15927, 21600, 22139, 22825, 27289, 27999, 28718, 29305, 31134, 31269, 32463, 33195, 36387, 36390, 38888])
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
H=set()
for i in broken: H|=(eqvars[i]&freeinp)
# for each handle, count how many equations OUTSIDE broken∪core it appears in (=ripple cost)
info=[]
for h in sorted(H):
    outside=byh[h]-broken-CORE
    info.append((len(outside), h, len(byh[h])))
info.sort()
print("handle | total_eqs | ripple(outside broken+core):")
for cost,h,tot in info[:40]:
    print(f"  x_{h}: total={tot}, ripple={cost}")
print(f"\nhandles with ripple==0 (private to broken set): {[h for c,h,t in info if c==0]}")
