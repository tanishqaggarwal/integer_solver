#!/usr/bin/env python3
"""Test wire inertness (methodology escape). Set all 220 wire members to sign*V for V in {1,2,p+1},
keep everything else at best, count fails. If inert -> few new fails -> the wire is a free parameter
and the core is solvable by product-slacks (x_5101*x_30317=-L1 with x_5101=V etc.)."""
import json, re, ast, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
A=load_atoms()
par=list(range(NVARS)); sgn=[1]*NVARS
def find2(x):
    s=1; r=x
    while par[r]!=r: s*=sgn[r]; r=par[r]
    return r,s
def union(a,b,rel):
    ra,sa=find2(a); rb,sb=find2(b)
    if ra==rb: return
    par[ra]=rb; sgn[ra]=rel*sb*sa
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==2 and pp.get((),0)==0:
        v1,v2=sorted(vs); c1=pp.get((v1,),0); c2=pp.get((v2,),0)
        qok=all(pp.get(k,0)==0 for k in pp if isinstance(k,tuple) and len(k)==2)
        if qok and c1!=0 and c2!=0 and abs(c1)==abs(c2):
            rel=(-c2)//c1
            if rel in (1,-1): union(v1,v2,rel)
r0,_=find2(26064)
wire={v:find2(v)[1] for v in range(NVARS) if find2(v)[0]==r0}
print(f"wire size {len(wire)}")
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
def count(setwire):
    val=[0]*NVARS
    for k,x in best.items():
        if k<NVARS: val[k]=x
    if setwire is not None:
        for v,s in wire.items(): val[v]=s*setwire
    ns={'__builtins__':{},'v':val}
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    return F,val
F0,_=count(None)
print(f"baseline (wire=p): {len(lines)-len(F0)} sat, {len(F0)} fail (core {len([i for i in F0 if i in CORE])})")
for V in [1,2]:
    F,val=count(V)
    nc=[i for i in F if i not in CORE]
    print(f"wire=sign*{V}: {len(lines)-len(F)} sat, {len(F)} fail (core {len([i for i in F if i in CORE])}, noncore {len(nc)}); noncore sample {sorted(nc)[:12]}")
