import heal_harness as H
import json
from collections import defaultdict
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val[:]
baseF=set(H.fails())
# free ancestors of the gate-side gap vars
for t in [2099,19964,24908,17601]:
    a=sorted(H.anc.get(t,set()))
    print(f"x_{t}: {len(a)} free ancestors: {a[:40]}")
print()
# For each free ancestor of x_2099 and x_19964, compute fan-out (satisfied eqs it's in) and current val
def report(name, t):
    print(f"=== free ancestors of x_{name} (shift changes {name}, hence G) ===")
    rows=[]
    for v in sorted(H.anc.get(t,set())):
        eqs=[i for i,vs in enumerate(H.eqvars) if v in vs]
        sat=[i for i in eqs if i not in baseF]
        rows.append((len(sat),v,H.val[v],len(eqs)))
    rows.sort()
    for nsat,v,val,neq in rows[:25]:
        print(f"  x_{v}: fanout(satisfied)={nsat}  neq={neq}  val={val}")
report('2099',2099)
report('19964',19964)
