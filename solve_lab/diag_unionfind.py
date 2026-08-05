#!/usr/bin/env python3
"""Signed union-find over ALL 2-term identity atoms (x_a - x_b and x_a + x_b). Find the class of
x_26064 (forced=p). Check if x_5101(→x_4376), x_32017(→x_16123), x_26789(→x_35148) are in it.
If SEPARATE, they are NOT forced to p — free to set to 1 (the methodology escape)."""
import json, re, ast, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
A=load_atoms()
# union-find with sign
par=list(range(NVARS)); sgn=[1]*NVARS
def find(x):
    root=x; s=1
    while par[root]!=root: s*=sgn[root]; root=par[root]
    # path compress
    cur=x; s2=1
    while par[cur]!=cur:
        nxt=par[cur]; ns=sgn[cur]; par[cur]=root; sgn[cur]=s2*s if False else None
        cur=nxt
    return root
# simpler find without compression bug:
par=list(range(NVARS)); sgn=[1]*NVARS
def find2(x):
    s=1; r=x
    while par[r]!=r: s*=sgn[r]; r=par[r]
    return r,s
def union(a,b,rel):  # x_a = rel * x_b  (rel=+1 or -1)
    ra,sa=find2(a); rb,sb=find2(b)
    if ra==rb: return
    # x_a=sa*root_a, x_b=sb*root_b, want x_a=rel*x_b -> sa*RA = rel*sb*RB -> RA = (rel*sb/sa) RB
    par[ra]=rb; sgn[ra]=rel*sb*sa  # since sa,sb in {+1,-1}, 1/sa=sa
# parse 2-term identity atoms from atoms
cnt=0
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==2 and pp.get((),0)==0:
        v1,v2=sorted(vs)
        c1=pp.get((v1,),0); c2=pp.get((v2,),0)
        # pure linear 2-term: c1*x1 + c2*x2 = 0, no quadratic
        qok=all(pp.get(k,0)==0 for k in pp if isinstance(k,tuple) and len(k)==2)
        if qok and c1!=0 and c2!=0 and abs(c1)==abs(c2):
            rel=-1 if (c1*c2>0) else 1  # c1 x1 + c2 x2=0 -> x1 = -c2/c1 x2
            rel = (-c2)//c1 if c1!=0 else 0
            if rel in (1,-1):
                union(v1,v2,rel); cnt+=1
print(f"identity unions: {cnt}")
r26064,_=find2(26064)
cls=[v for v in range(NVARS) if find2(v)[0]==r26064]
print(f"x_26064 class size: {len(cls)}")
for v in [5101,4376,32017,16123,26789,35148]:
    rv,sv=find2(v)
    print(f"  x_{v}: root={rv}, in x_26064 class: {rv==r26064}, sign={sv}")
