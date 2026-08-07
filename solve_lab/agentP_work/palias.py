#!/usr/bin/env python3
"""Agent P: structural constant classes (alias closure), independent of any assignment."""
import pickle,sys,json
from collections import defaultdict,Counter
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
NV=38748
par=list(range(NV))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def uni(a,b):
    a,b=find(a),find(b)
    if a!=b: par[a]=b
pins={}
for ap in AP:
    ks=list(ap.items())
    if len(ks)==2:
        (m1,c1),(m2,c2)=sorted(ks,key=lambda z:len(z[0]))
        if len(m1)==1 and len(m2)==1 and c1==-c2: uni(m1[0],m2[0])
        elif m1==() and len(m2)==1 and abs(c2)==1:
            pins.setdefault(m2[0],set()).add(-c1//c2)
    if len(ks)==1 and list(ap)[0]!=() and len(list(ap)[0])==1:
        pins.setdefault(list(ap)[0][0],set()).add(0)
cls=defaultdict(set)
for x in range(NV): cls[find(x)].add(x)
# propagate pins through classes
clspin=defaultdict(set)
for x,vs in pins.items(): clspin[find(x)] |= vs
big=[(len(cls[r]),r,clspin.get(r,set())) for r in cls if clspin.get(r)]
big.sort(reverse=True)
print("constant classes (size, pin values):")
for n,r,v in big[:12]:
    vv={ (str(x)[:26]+'..%dbit'%x.bit_length()) if abs(x)>10**9 else x for x in v}
    print(f"  size={n} pins={vv}")
out={}
for n,r,v in big:
    for x in cls[r]: out[x]=sorted(v)[0] if len(v)==1 else None
pickle.dump({'const':out,'par':[find(x) for x in range(NV)],'clspin':{k:sorted(v) for k,v in clspin.items()}},open(W+'alias.pkl','wb'))
print("total vars pinned to a constant:",sum(1 for x in out if out[x] is not None))
