#!/usr/bin/env python3
"""Extract the 13 unpacking + codeword eqs as polynomials in the 9 bits at wire=1 (else agentA).
Check separability (cross-terms) and structure."""
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
UNP=[8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
CODE=[8135,25994,19035,19805]  # codeword checks
ns={'v':v,'__builtins__':{}}
def ev(i): return eval(eqcode[i],ns)
# save bit values, extract
saved={b:v[b] for b in BITS}
def setbits(vals):
    for b in BITS: v[b]=vals.get(b,0)
# constant (all bits 0)
def poly_of(i):
    setbits({})
    c0=ev(i)
    lin={}; quad={}
    for b in BITS:
        if b not in eqvars[i]: continue
        setbits({b:1}); f1=ev(i)
        setbits({b:2}); f2=ev(i)
        setbits({})
        # f(t)=c0+lin*t+quad*t^2 ; f1=c0+lin+quad; f2=c0+2lin+4quad
        quad_b=(f2-2*f1+c0)//2
        lin_b=f1-c0-quad_b
        if lin_b or quad_b: lin[b]=lin_b; quad[b]=quad_b
    # cross terms
    cross={}
    bs=[b for b in BITS if b in eqvars[i]]
    for a_i in range(len(bs)):
        for b_i in range(a_i+1,len(bs)):
            ba,bb=bs[a_i],bs[b_i]
            setbits({ba:1,bb:1}); fab=ev(i); setbits({})
            cr=fab-c0-lin.get(ba,0)-quad.get(ba,0)-lin.get(bb,0)-quad.get(bb,0)
            if cr: cross[(ba,bb)]=cr
    return c0,lin,quad,cross
print("=== 13 unpacking eqs at wire=1 (poly in 9 bits) ===")
for i in UNP:
    c0,lin,quad,cross=poly_of(i)
    print(f"eq {i}: const%p={c0%p}, quad coeffs={ {b:quad[b] for b in quad} }, cross={len(cross)} terms")
print("\n=== codeword eqs ===")
for i in CODE:
    if i>=len(lines): continue
    c0,lin,quad,cross=poly_of(i)
    print(f"eq {i}: fails-now={ev(i)!=0}, const%p={c0%p}, #quad={len(quad)}, #lin={len(lin)}, cross={len(cross)}")
setbits(saved)
