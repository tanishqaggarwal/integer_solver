#!/usr/bin/env python3
"""Extract ALL distinct gadget conditions from the 27 (0,0) failures. Group them."""
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
# collect distinct nonzero gadget sources
gadgets=defaultdict(list)
for i in F:
    node=inner(lines[i].rsplit('=',1)[0])
    for s,t in flat(node):
        if s*evn(t)!=0:
            src=ast.unparse(t)
            gadgets[src].append(i)
print(f"{len(F)} fails, {len(gadgets)} distinct nonzero gadgets:")
for src,eqs in sorted(gadgets.items(), key=lambda x:-len(x[1])):
    v=eval(compile(VAR.sub(r'v[\1]',src),'<e>','eval'),ns)
    fv=[int(m) for m in VAR.findall(src) if int(m) in freeinp]
    print(f"  [{v.bit_length()}b in {len(eqs)} eqs] {src[:70]}  free={fv}")
