#!/usr/bin/env python3
"""Deep structural analysis of the 20 quadratic core equations. For each, find the inner squared
expression E, decompose E = sum of gadget terms, and identify which free inputs control it and
whether those free inputs are shared with (coupled to) the linear equations."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)

# load best partial solution
sol={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
def V(v): return sol.get(v,0)

lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
CORE=[2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]

VAR=re.compile(r'x_(\d+)')
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

# expand a gate output symbolically down to free inputs, tracking the actual definition tree
def expand_def(t, depth=0):
    """Return string description of what gate t computes, one level."""
    if t in freeinp: return f"free x_{t}(={V(t)})"
    rhs,vids=gate_defs[t]
    return rhs

for i in CORE:
    lhs=lines[i].rsplit('=',1)[0]
    node=inner(lhs)
    src=ast.unparse(node)
    vs=set(int(m) for m in VAR.findall(src))
    print(f"\n=== eq {i}: E = {src[:120]}")
    print(f"    vars in E: {sorted(vs)}")
    for v in sorted(vs):
        tag = "FREE" if v in freeinp else "gate"
        d = "" if v in freeinp else f" := {gate_defs[v][0][:80]}"
        print(f"      x_{v} [{tag}] = {V(v)}{d}")
