#!/usr/bin/env python3
"""For each failing equation, break into gadgets (out - expr) and find which gadget is nonzero
under the best partial. Pinpoints the exact failing gate + whether it's a multi-def join point."""
import json, re, ast
from collections import defaultdict
from propagate import NVARS
cand={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39007.json')).items()}
val=[0]*NVARS
for k,x in cand.items():
    if k<NVARS: val[k]=x
mt=json.load(open('atoms/multi_targets.json'))  # dict var->ndefs
multidef=set(int(k) for k in mt)
# map each var to its gate defs
gate_defs=defaultdict(list)
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gate_defs[d['t']].append((d['rhs'],d['eq']))
fails=json.load(open('fails_39007.json'))
lines=[L for L in open('EQUATIONS.txt'.replace('EQUATIONS','../EQUATIONS')).read().split('\n') if L.strip()]

def ev(node):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -ev(node.operand)
    if isinstance(node,ast.BinOp):
        a=ev(node.left); b=ev(node.right)
        if isinstance(node.op,ast.Add): return a+b
        if isinstance(node.op,ast.Sub): return a-b
        if isinstance(node.op,ast.Mult): return a*b
    raise ValueError

# strip square/const to get the inner E
def inner(lhs):
    node=ast.parse(lhs,mode='eval').body
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,ast.Constant) or (isinstance(a,ast.UnaryOp) and isinstance(a.operand,ast.Constant))
        cb=isinstance(b,ast.Constant) or (isinstance(b,ast.UnaryOp) and isinstance(b.operand,ast.Constant))
        if ca and not cb: node=b
        elif cb and not ca: node=a
        elif (not ca) and (not cb) and ast.unparse(a)==ast.unparse(b): node=a
        else: break
    return node

# flatten E into additive gadget terms (coef * gadget)
def flatten(node, sign=1, out=None):
    if out is None: out=[]
    if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub)):
        flatten(node.left, sign, out)
        flatten(node.right, sign*(1 if isinstance(node.op,ast.Add) else -1), out)
    else:
        out.append((sign,node))
    return out

print("Localizing failing gadgets (culprit gate + join-point status):\n")
culprits=defaultdict(int)
for i in fails:
    E=inner(lines[i].rsplit("=",1)[0])
    terms=flatten(E)
    bad=[]
    for sgn,t in terms:
        v=ev(t)
        if v!=0:
            # identify the output var of this gadget (lone name on the left of a subtraction)
            src=ast.unparse(t)
            m=re.match(r'\(?x_(\d+)\)?\s*-', src.replace(' ',''))
            outv=int(m.group(1)) if m else None
            bad.append((outv, v.bit_length(), src[:45]))
    tag=[]
    for outv,bl,src in bad:
        j = 'JOIN' if outv in multidef else ('input' if outv is not None and outv not in gate_defs else 'gate')
        tag.append(f"x_{outv}[{j},2^{bl}]")
        if outv is not None: culprits[outv]+=1
    print(f"eq[{i}]: {len(bad)} bad gadget(s): {', '.join(tag)}")
print("\nCulprit output vars (appear as failing gadget across eqs):")
for v,n in sorted(culprits.items(), key=lambda x:-x[1]):
    j='JOIN' if v in multidef else ('input' if v not in gate_defs else 'gate')
    print(f"  x_{v}: {n} eqs [{j}], defs={len(gate_defs.get(v,[]))}")
