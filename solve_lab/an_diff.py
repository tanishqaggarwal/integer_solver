#!/usr/bin/env python3
import re,json
p=2**256-2**32-977; NVARS=38748
VAR_RE=re.compile(r'x_(\d+)')
FAILS=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
def loadd(path):
    d=json.load(open(path)); out={}
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); out[idx]=int(val)
    return out
a=loadd('best/new_instance_partial_39013.json')
b=loadd('best_agentA_39022.json')
keys=set(a)|set(b)
changed=[k for k in keys if a.get(k,0)!=b.get(k,0)]
print(f"changed vars (39013 -> agentA 39022): {len(changed)}")
print("changed:", sorted(changed)[:60])
# which of the 11's variables were changed?
codes={}; varsets={}
with open('../EQUATIONS.txt') as f:
    for i,line in enumerate(f):
        line=line.strip()
        if i in FAILS:
            lhs=line.rsplit('=',1)[0]
            varsets[i]=set(int(m) for m in VAR_RE.findall(lhs))
allv=set()
for i in FAILS: allv|=varsets[i]
changed_in_11=[k for k in changed if k in allv]
print(f"\nchanged vars that appear in the 11: {sorted(changed_in_11)}")
# per-equation: which changed vars does each of the 11 contain
for i in FAILS:
    ch=sorted(varsets[i]&set(changed))
    print(f"  eq {i}: changed vars present = {ch}")
