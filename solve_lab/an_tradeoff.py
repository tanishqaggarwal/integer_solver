#!/usr/bin/env python3
import re,json
p=2**256-2**32-977; NVARS=38748
VAR_RE=re.compile(r'x_(\d+)')
FAILS=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
# core square eqs: find the 20. They are E^2 form. Let me detect by checking known.
# Load all eqs compiled
codes=[]; 
with open('../EQUATIONS.txt') as f:
    for line in f:
        line=line.strip()
        if not line: codes.append(None); continue
        lhs=line.rsplit('=',1)[0]
        codes.append(compile(VAR_RE.sub(r'v[\1]',lhs),'<eq>','eval'))
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
def check(v,idxs):
    ns={'v':v,'__builtins__':{}}
    return [i for i in idxs if eval(codes[i],ns)!=0]
for name in ['best/new_instance_partial_39013.json','best/new_instance_partial_39018.json','best/new_instance_partial_39021.json','best_agentA_39022.json']:
    try: v=loadv(name)
    except Exception as e: print(name,"ERR",e); continue
    f11=check(v,FAILS)
    # total fails
    allf=check(v,range(len(codes)))
    print(f"{name}: total_fail={len(allf)}  of-the-11-fail={len(f11)} -> {f11}")
