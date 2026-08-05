#!/usr/bin/env python3
"""Analyze the 11 residual failing equations: mod-p residuals at agentA solution,
shared-variable structure, and whether they reduce to few scalar conditions."""
import re, json, time
p = 2**256 - 2**32 - 977
NVARS=38748
VAR_RE=re.compile(r'x_(\d+)')
FAILS=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]

# compile equations
codes={}; varsets={}
with open('../EQUATIONS.txt') as f:
    for i,line in enumerate(f):
        line=line.strip()
        if not line: continue
        if i in FAILS:
            lhs=line.rsplit('=',1)[0]
            ids=tuple(int(m) for m in VAR_RE.findall(lhs))
            codes[i]=compile(VAR_RE.sub(r'v[\1]',lhs),'<eq>','eval')
            varsets[i]=set(ids)

# load agentA solution
d=json.load(open('best_agentA_39022.json'))
v=[0]*NVARS
for k,val in d.items():
    idx=int(k[2:]) if k.startswith('x_') else int(k)
    v[idx]=int(val)

ns={'v':v,'__builtins__':{}}
print("=== residuals mod p at agentA solution ===")
res={}
for i in FAILS:
    r=eval(codes[i],ns)
    res[i]=r%p
    print(f"eq {i}: resid%p = {r%p}   (exact int nonzero: {r!=0})")

# variable sharing among the 11
from collections import Counter
c=Counter()
for i in FAILS:
    for x in varsets[i]: c[x]+=1
shared=[(x,n) for x,n in c.items() if n>=2]
shared.sort(key=lambda t:-t[1])
print(f"\n=== variables shared by >=2 of the 11 (total distinct vars {len(c)}) ===")
print("top shared:", shared[:30])

# union of all variables
allv=set()
for i in FAILS: allv|=varsets[i]
print(f"\ntotal distinct variables across 11 eqs: {len(allv)}")
