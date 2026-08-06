#!/usr/bin/env python3
"""For each of the 13 eqs that break when wire=1, flat-decompose the LHS and find the terms whose
value changes p->1. Identify the wire-partner in each (the handle to compensate)."""
import json, re, ast
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
r0,_=find2(26064); wire={v:find2(v)[1] for v in range(NVARS) if find2(v)[0]==r0}
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
def mkval(W):
    val=[0]*NVARS
    for k,x in best.items():
        if k<NVARS: val[k]=x
    for v,s in wire.items(): val[v]=s*W
    return val
def flat(node,s=1,o=None):
    if o is None:o=[]
    if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub)):
        flat(node.left,s,o); flat(node.right,s*(1 if isinstance(node.op,ast.Add) else -1),o)
    else: o.append((s,node))
    return o
def evn(node,val):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -evn(node.operand,val)
    a=evn(node.left,val); b=evn(node.right,val)
    return a+b if isinstance(node.op,ast.Add) else a-b if isinstance(node.op,ast.Sub) else a*b
b13=[8429, 11166, 11915, 12594, 23869, 25313, 26785, 31400, 32300, 36106, 36767, 37257]
vp=mkval(p); v1=mkval(1)
for i in b13[:6]:
    node=ast.parse(lines[i].rsplit('=',1)[0],mode='eval').body
    # descend to the innermost non-const-mult
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,(ast.Constant,)); cb=isinstance(b,(ast.Constant,))
        if ca: node=b
        elif cb: node=a
        else: break
    terms=flat(node)
    print(f"\neq {i}: changed terms (p->1):")
    for s,t in terms:
        dp=evn(t,vp); d1=evn(t,v1)
        if dp!=d1:
            src=ast.unparse(t)
            wv=[int(m) for m in VAR.findall(src) if int(m) in wire]
            nz=[int(m) for m in VAR.findall(src) if best.get(int(m),0)!=0 and int(m) not in wire]
            print(f"   Δ [{src[:55]}] wire={wv} nonzero_nonwire={nz}")
