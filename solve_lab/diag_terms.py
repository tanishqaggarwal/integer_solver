#!/usr/bin/env python3
"""For each core eq, flatten inner E into additive terms; evaluate each term at the best partial.
Show nonzero terms and the free inputs controlling them."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
sol={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
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
CORE=[2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]
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
def freevars_of(node):
    return set(int(m) for m in VAR.findall(ast.unparse(node)) if int(m) in freeinp)
for i in CORE:
    node=inner(lines[i].rsplit('=',1)[0])
    E=evn(node)
    terms=flat(node)
    nz=[(s,t) for s,t in terms if s*evn(t)!=0]
    print(f"\n=== eq {i}: E={E} ({'ZERO' if E==0 else 'bits='+str(E.bit_length())}), {len(terms)} terms, {len(nz)} nonzero")
    for s,t in nz:
        tv=evn(t); src=ast.unparse(t)
        fvs=freevars_of(t)
        print(f"    {'+' if s>0 else '-'} [{tv.bit_length()}b] {src[:75]}  free={sorted(fvs)}")
