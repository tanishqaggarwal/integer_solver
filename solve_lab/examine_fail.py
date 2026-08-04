#!/usr/bin/env python3
"""Ground-truth: which RAW equations fail under best partial, and their structure."""
import json, re, ast
from propagate import NVARS

# load best partial
cand = {int(k[2:]): v for k,v in json.load(open('best/new_instance_partial_39007.json')).items()}
val = [0]*NVARS
for v,x in cand.items(): val[v]=x

lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
print(f"{len(lines)} equations, NVARS={NVARS}")

VAR=re.compile(r'x_(\d+)')
def ev_ast(node):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -ev_ast(node.operand)
    if isinstance(node,ast.BinOp):
        a=ev_ast(node.left); b=ev_ast(node.right)
        if isinstance(node.op,ast.Add): return a+b
        if isinstance(node.op,ast.Sub): return a-b
        if isinstance(node.op,ast.Mult): return a*b
    raise ValueError

def is_square(lhs):
    """Detect top-level perfect square C*(E)*(E) or (E)*(E)."""
    node=ast.parse(lhs,mode='eval').body
    # strip leading constant multipliers
    def strip(n):
        while isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult):
            a,b=n.left,n.right
            ca=isinstance(a,ast.Constant)
            cb=isinstance(b,ast.Constant)
            if ca and not cb: n=b
            elif cb and not ca: n=a
            else: break
        return n
    n=strip(node)
    if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult):
        if ast.unparse(n.left)==ast.unparse(n.right):
            return n.left
    return None

fails=[]
for i,L in enumerate(lines):
    lhs=L.rsplit('=',1)[0]
    node=ast.parse(lhs,mode='eval').body
    r=ev_ast(node)
    if r!=0:
        ids=sorted(set(int(m) for m in VAR.findall(L)))
        sq=is_square(lhs)
        sqroot=None
        if sq is not None:
            sqroot=ev_ast(sq)
        fails.append((i,r,len(L),len(ids),sq is not None,sqroot))

print(f"\n{len(fails)} FAILING equations:")
for i,r,ll,nids,issq,sqr in fails:
    tag="SQUARE" if issq else "linear/other"
    extra=f" root={sqr}" if issq else ""
    print(f"  eq[{i}]: |resid|~2^{r.bit_length()-1 if r else 0} nvars={nids} len={ll} {tag}{extra}")
json.dump([f[0] for f in fails], open('fails_39007.json','w'))
