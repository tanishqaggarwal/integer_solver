#!/usr/bin/env python3
"""Decompose the (0,0)-quadrant failing equations at the q00_partial solution. See if they reduce
to a simpler condition than (1,1)'s S,T=0 (e.g. just quotient handles or a smaller monster set)."""
import json, re, ast, sys
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
sol={int(k[2:]):v for k,v in json.load(open('q00_partial.json')).items()}
val=[0]*NVARS
for v,x in sol.items():
    if v<NVARS: val[v]=x
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
def evn(node):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -evn(node.operand)
    a=evn(node.left); b=evn(node.right)
    return a+b if isinstance(node.op,ast.Add) else a-b if isinstance(node.op,ast.Sub) else a*b
def inner(lhs):
    node=ast.parse(lhs,mode='eval').body
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,ast.Constant) or (isinstance(a,ast.UnaryOp) and isinstance(a.operand,ast.Constant))
        cb=isinstance(b,ast.Constant) or (isinstance(b,ast.UnaryOp) and isinstance(b.operand,ast.Constant))
        if ca and not cb: node=b
        elif cb and not ca: node=a
        elif ast.unparse(a)==ast.unparse(b): node=a
        else: break
    return node
def flat(node,s=1,o=None):
    if o is None:o=[]
    if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub)):
        flat(node.left,s,o); flat(node.right,s*(1 if isinstance(node.op,ast.Add) else -1),o)
    else: o.append((s,node))
    return o
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{'v':val}}; ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"q00 partial fails: {len(F)}: {F}")
print(f"x_15298={val[15298]}")
# for each failing eq, show nonzero flat terms
for i in F[:8]:
    node=inner(lines[i].rsplit('=',1)[0]); E=evn(node)
    terms=flat(node); nz=[(s,t) for s,t in terms if s*evn(t)!=0]
    print(f"\neq {i}: E bits={E.bit_length() if E else 0}, {len(nz)} nonzero terms")
    for s,t in nz[:6]:
        src=ast.unparse(t); tv=evn(t)
        fv=[int(m) for m in VAR.findall(src) if int(m) in freeinp]
        print(f"   [{tv.bit_length()}b] {src[:60]} free={fv}")
