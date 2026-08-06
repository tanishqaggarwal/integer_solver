import json, re
from collections import defaultdict
import heal_harness as H
p=2**256-2**32-977
with open('atoms/poly_atoms.jsonl') as f:
    lines=f.readlines()
# parse pins: sel*(target - CONST) - coef*handle
pinrec=[]  # (atomidx, sel, target, const, coef, handle)
pat=re.compile(r'x_(\d+) \* \(x_(\d+) - (\d+)\) - (?:(\d+) \* )?x_(\d+)')
for i,l in enumerate(lines):
    a=json.loads(l)
    m=pat.match(a['repr'])
    if m:
        sel=int(m.group(1)); tgt=int(m.group(2)); const=int(m.group(3))
        coef=int(m.group(4)) if m.group(4) else 1; handle=int(m.group(5))
        pinrec.append((i,sel,tgt,const,coef,handle))
print(f"parsed {len(pinrec)} pins")
json.dump(pinrec, open('pinrec.json','w'))
# group by selector
bysel=defaultdict(list)
for rec in pinrec: bysel[rec[1]].append(rec)
# which bits feed the core? core control gates x_3558=x_24908-x_16742, x_29322=x_14853-x_12186
coretargets={24908,16742,14853,12186, 6418,31861,14865,12553, 2099,19964,7068,4432}
print("\ncore-relevant pin targets:")
for rec in pinrec:
    if rec[2] in coretargets:
        print(f"  atom{rec[0]}: sel=x_{rec[1]} loads x_{rec[2]}={rec[3]%p} (mod p), handle=x_{rec[5]}")
