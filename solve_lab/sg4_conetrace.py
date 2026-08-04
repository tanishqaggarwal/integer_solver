import heal_harness as H
import json, ast, re
from collections import defaultdict
p=H.p
# gate map: target -> (rhs_ast, vids)
gmap={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line)
        gmap[d['t']]=d['rhs']
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val
def is_wire(v): return V[v]!=0 and V[v]%p==0
# recursively expand target through +/- into leaves; classify leaves
slacks=[]  # (product var, factorA, factorB)
visited=set()
def expand(t, depth=0):
    if t in visited or depth>40: return
    visited.add(t)
    if t not in gmap:  # free input
        return
    rhs=gmap[t]
    node=ast.parse(rhs,mode='eval').body
    # walk additive structure
    def walk(n):
        if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub)):
            walk(n.left); walk(n.right); return
        if isinstance(n,ast.UnaryOp):
            walk(n.operand); return
        # leaf term: could be Name, Mult, Const
        if isinstance(n,ast.Name):
            expand(int(n.id[2:]),depth+1); return
        if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult):
            # product term
            vs=[int(m.group(1)) for m in re.finditer(r'x_(\d+)',ast.unparse(n))]
            # is it a genuine product slack (2 vars, not const*var)?
            names=[x for x in [n.left,n.right]]
            varnames=re.findall(r'x_(\d+)',ast.unparse(n))
            if len(set(varnames))>=2:
                a,b=int(varnames[0]),int(varnames[-1])
                slacks.append((t,a,b))
            for vv in vs: expand(vv,depth+1)
            return
    walk(node)
for target in [19964, 2099]:
    slacks.clear(); visited.clear()
    expand(target)
    print(f"\n===== product-slacks in additive cone of x_{target} ({len(slacks)}) =====")
    seen=set()
    for host,a,b in slacks:
        key=(a,b)
        if key in seen: continue
        seen.add(key)
        # classify: which factor is free & currently 0 (activatable), which is the multiplier
        def cls(v):
            return f"x_{v}[free={v in H.freeinp},wire={is_wire(v)},val={V[v] if abs(V[v])<1e10 else ('0' if V[v]==0 else 'BIG')}]"
        # determine gating of the multiplier factor: trace if it depends on a selector bit
        print(f"  slack in x_{host}: {cls(a)} * {cls(b)}")
