#!/usr/bin/env python3
"""Find all multiply-gate wire-products (wire_member * partner) with nonzero partner. Report each
partner, its value, and how many equations it appears in (privacy). This is the set to scale by p
when freeing the wire (except the 3 core partners x_30317,x_5146,x_2936)."""
import json, re
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
r0,_=find2(26064); wire=set(v for v in range(NVARS) if find2(v)[0]==r0)
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
def V(v): return best.get(v,0)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
occ=defaultdict(int)
for L in lines:
    for m in set(int(x) for x in VAR.findall(L)): occ[m]+=1
# find multiply gates with a wire factor
partners=[]
for t,rhs,vids in gates:
    m=re.fullmatch(r'x_(\d+) \* x_(\d+)', rhs)
    if m:
        a,b=int(m.group(1)),int(m.group(2))
        wf = a in wire, b in wire
        if wf[0]^wf[1]:  # exactly one is wire
            partner = b if wf[0] else a
            if V(partner)!=0:
                partners.append((partner, t, V(partner).bit_length(), occ[partner]))
core={30317,5146,2936}
print(f"wire size {len(wire)}; nonzero wire-product partners: {len(partners)}")
print("partner -> (gate, partner_bits, #eqs_appearing):")
for pt,t,bits,oc in sorted(partners, key=lambda x:-x[3]):
    tag=" [CORE]" if pt in core else ""
    print(f"  x_{pt}: gate x_{t}, {bits} bits, in {oc} eqs{tag}")
