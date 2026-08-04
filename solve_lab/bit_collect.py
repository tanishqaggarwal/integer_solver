#!/usr/bin/env python3
"""Collect ALL equations involving the 13 bits; extract bit-polynomials at wire=1 (else agentA).
Report the full constraint set size and structure."""
import json,pickle,re
from propagate import NVARS
p=2**256-2**32-977
D=pickle.load(open('wire_data.pkl','rb')); wire=D['wire']
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
v=vA[:]
for y,s in wire.items(): v[y]=s*1
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[frozenset(int(m) for m in VAR.findall(L)) for L in lines]
BITS=[31342,32058,22473,17389,1488,28827,37384,11094,18211,875,29159,37076,14048]
BITSET=set(BITS)
# all eqs containing any bit
BE=[i for i in range(len(lines)) if eqvars[i]&BITSET]
print(f"equations containing >=1 of the 13 bits: {len(BE)}")
ns={'v':v,'__builtins__':{}}
def ev(i): return eval(eqcode[i],ns)
# how many currently fail at wire=1 (bits=agentA)?
failing=[i for i in BE if ev(i)!=0]
print(f"  of those, failing at wire=1(bits=agentA): {len(failing)} -> {failing}")
# how many involve ONLY bits + wire (pure bit eqs given wire fixed)?
def setbits(vals):
    for b in BITS: v[b]=vals.get(b,0)
saved={b:v[b] for b in BITS}
# For each bit-eq, extract poly in bits; check separability and coeff magnitudes
polys={}
huge=[]
for i in BE:
    setbits({}); c0=ev(i)
    lin={};quad={}
    bs=[b for b in BITS if b in eqvars[i]]
    for b in bs:
        setbits({b:1});f1=ev(i);setbits({b:2});f2=ev(i);setbits({})
        q=(f2-2*f1+c0)//2; l=f1-c0-q
        quad[b]=q; lin[b]=l
    cross={}
    for a_ in range(len(bs)):
        for b_ in range(a_+1,len(bs)):
            ba,bb=bs[a_],bs[b_]
            setbits({ba:1,bb:1});fab=ev(i);setbits({})
            cr=fab-c0-lin[ba]-quad[ba]-lin[bb]-quad[bb]
            if cr: cross[(ba,bb)]=cr
    polys[i]={'c0':c0,'lin':lin,'quad':quad,'cross':cross,'bits':bs}
    if any(abs(q)>10**30 for q in quad.values()): huge.append(i)
setbits(saved)
print(f"  eqs with huge (>1e30) quad coeffs: {len(huge)} -> {huge}")
# separability
nsep=[i for i in BE if polys[i]['cross']]
print(f"  eqs with cross-terms (non-separable): {len(nsep)} -> {nsep}")
pickle.dump({'polys':polys,'BE':BE,'BITS':BITS,'failing':failing}, open('bitpolys.pkl','wb'))
print("saved bitpolys.pkl")
