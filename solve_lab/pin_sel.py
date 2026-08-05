import json, re
from collections import defaultdict, Counter
import heal_harness as H
p=2**256-2**32-977
with open('atoms/poly_atoms.jsonl') as f:
    lines=f.readlines()
VAR=re.compile(r'x_(\d+)')
pins=[]
for i,l in enumerate(lines):
    a=json.loads(l); poly=a['poly']
    bigcoef=[t for t in poly if t[0] and abs(t[1])>p//4]
    if bigcoef and len(poly)<=5:
        # parse: sel*(x_i - CONST) - coef*handle. selector = the var multiplied with x_i in a degree-2 term
        # find degree-2 term [selvar, ivar] with the big const being sel*(-CONST)
        deg2=[t for t in poly if len(t[0])==2]
        pins.append((i,a,deg2))
# selectors: for each pin, the selector is the var that appears in BOTH the deg2 product and alone with big const
sels=[]
targets=[]
for i,a,deg2 in pins:
    r=a['repr']
    m=re.match(r'x_(\d+) \* \(x_(\d+) -', r)
    if m:
        sels.append(int(m.group(1))); targets.append(int(m.group(2)))
print(f"pins parsed: {len(sels)}")
sc=Counter(sels)
print(f"distinct selectors: {len(sc)}")
print(f"selector frequency (top): {sc.most_common(10)}")
# are selectors free or gates?
freesel=[s for s in sc if s in H.freeinp]
gatesel=[s for s in sc if s not in H.freeinp]
print(f"selectors that are FREE inputs: {len(freesel)}: {sorted(freesel)[:20]}")
print(f"selectors that are GATES: {len(gatesel)}: {sorted(gatesel)[:20]}")
# current values of selectors (at sy config)
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
selvals=Counter(H.val[s] for s in sc)
print(f"selector values at sy config: {dict(list(selvals.items())[:8])}")
# targets that are free
freetgt=[t for t in targets if t in H.freeinp]
print(f"pin targets (message vars): {len(set(targets))} distinct, {len(set(freetgt))} free")
