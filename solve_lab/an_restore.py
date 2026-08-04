#!/usr/bin/env python3
import re,json,time
p=2**256-2**32-977; NVARS=38748
VAR_RE=re.compile(r'x_(\d+)')
FAILS=set([2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125])
t0=time.time()
codes=[]
with open('../EQUATIONS.txt') as f:
    for line in f:
        line=line.strip()
        if not line: codes.append(None); continue
        codes.append(compile(VAR_RE.sub(r'v[\1]',line.rsplit('=',1)[0]),'<eq>','eval'))
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
v013=loadv('best/new_instance_partial_39013.json')
def allfails(v):
    ns={'v':v,'__builtins__':{}}
    return set(i for i,c in enumerate(codes) if c is not None and eval(c,ns)!=0)
baseF=allfails(vA)
print(f"agentA base fails: {len(baseF)} loaded in {time.time()-t0:.1f}s")
# restore x_4432 and x_7068
v2=vA[:]
v2[4432]=v013[4432]; v2[7068]=v013[7068]
F2=allfails(v2)
print(f"after restoring x_4432,x_7068 to 39013 values: {len(F2)} fails")
print(f"  fixed (were failing, now ok): {sorted(baseF-F2)}")
print(f"  newly broken: {sorted(F2-baseF)}")
