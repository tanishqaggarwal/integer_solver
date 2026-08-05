import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('g1g2_closed.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
print(f"start: {len(F0)} fails")
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
# atom7450: x_2964 - x_26756 - x_13859*x_19569 = 0 -> x_19569 = (x_2964-x_26756)/x_13859
num1=V[2964]-V[26756]; den1=V[13859]
print(f"x_13859={V[13859]}, (x_2964-x_26756)%x_13859 == 0? {num1%den1==0 if den1 else 'div0'}")
# atom7452: 9367949*(x_24548-x_25442) - x_15616*x_11052 = 0 -> x_11052 = 9367949*(x_24548-x_25442)/x_15616
num2=9367949*(V[24548]-V[25442]); den2=V[15616]
print(f"x_15616={V[15616]}, num2%x_15616 == 0? {num2%den2==0 if den2 else 'div0'}")
# check x_19569, x_11052 are the free partners (private)
print(f"x_19569 {'FREE' if 19569 in H.freeinp else 'gate'}, x_11052 {'FREE' if 11052 in H.freeinp else 'gate'}")
print(f"x_13859 {'FREE' if 13859 in H.freeinp else 'gate'}, x_15616 {'FREE' if 15616 in H.freeinp else 'gate'}")
ok=True
if den1 and num1%den1==0: V[19569]=num1//den1
else: ok=False; print("  DIV1 fail")
if den2 and num2%den2==0: V[11052]=num2//den2
else: ok=False; print("  DIV2 fail")
if ok:
    fwd_from([19569,11052])
    print(f"atom7450 = {V[2964]-V[26756]-V[579]}")
    print(f"atom7452 = {9367949*(V[24548]-V[25442])-V[7927]}")
    F=set(H.fails())
    print(f"\n*** FAILS: {len(F)}: {sorted(F)} ***")
    print(f"fixed: {sorted(F0-F)}")
    print(f"broken: {sorted(F-F0)}")
    if len(F)==0:
        json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('SOLVED_full.json','w'))
        print("=== SOLVED_full.json SAVED ===")
    elif len(F)<11:
        json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('heal_slacks_partial.json','w'))
