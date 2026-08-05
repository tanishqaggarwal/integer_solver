#!/usr/bin/env python3
"""Decompose the 12 new failing equations in agentA's 39021. Find their nonzero gadget terms and
whether they reduce to a small set of conditions (like the original core reduced to M1,M2,M3)."""
import json, re, ast
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
s={int(k[2:]):v for k,v in json.load(open('best_agentA_39021.json')).items()}
val=[0]*NVARS
for k,v in s.items():
    if k<NVARS: val[k]=v
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
F12=[4833, 4944, 5348, 9344, 10406, 11574, 12321, 19708, 20927, 21972, 27514, 38014]
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
def flat(node,sg=1,o=None):
    if o is None:o=[]
    if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub)):
        flat(node.left,sg,o); flat(node.right,sg*(1 if isinstance(node.op,ast.Add) else -1),o)
    else: o.append((sg,node))
    return o
gadgets=defaultdict(list)
for i in F12:
    node=inner(lines[i].rsplit('=',1)[0])
    for sg,t in flat(node):
        if sg*evn(t)!=0:
            src=ast.unparse(t)
            gadgets[src].append(i)
print(f"12 fails, {len(gadgets)} distinct nonzero gadgets:")
for src,eqs in sorted(gadgets.items(),key=lambda x:-len(x[1]))[:20]:
    v=evn(ast.parse(src,mode='eval').body)
    fv=[int(m) for m in VAR.findall(src) if int(m) in freeinp]
    print(f"  [{v.bit_length()}b x{len(eqs)}] {src[:65]} free={fv}")
