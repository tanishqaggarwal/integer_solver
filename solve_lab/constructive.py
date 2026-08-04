#!/usr/bin/env python3
"""Constructive load-absorption for the FEASIBLE (1,1) config. Acyclic circuit => set each free
input to exactly zero its defining gadget, in forward-eval order. Repeat to convergence."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
setfree=set(override)  # free inputs we've deliberately set (don't overwrite as absorbers unless load)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp or pinned[v]: ready[v]=True
gu=[0]*len(gates); using=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates):
    u=0
    for v in vids:
        if not ready[v]: u+=1
        using[v].append(gi)
    gu[gi]=u
definer={}; order=[]
q=deque(gi for gi in range(len(gates)) if gu[gi]==0)
while q:
    gi=q.popleft(); t,rhs,vids=gates[gi]
    if ready[t]: continue
    definer[t]=gi; order.append(t); ready[t]=True
    for gj in using[t]:
        gu[gj]-=1
        if gu[gj]==0: q.append(gj)
VAR=re.compile(r'x_(\d+)')
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
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
    else: o.append(node)
    return o
astcache={}
def rootast(i):
    if i not in astcache: astcache[i]=inner(lines[i].rsplit('=',1)[0])
    return astcache[i]
# Try to set a free input to make gadget-term t == 0. t is coef*(subexpr). Find a free input v in
# subexpr that appears linearly with a computable coefficient; set val[v] so subexpr's value 0.
def try_absorb(t):
    # descend into t to find the "gadget core" (strip constant multipliers and set activators)
    node=t
    # find free inputs in node
    names=[int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node))]
    freevars=[v for v in names if v in freeinp]
    # evaluate d(node)/d(v) by finite diff; pick v with nonzero deriv and |current node| divisible
    cur=evn(node)
    if cur==0: return False
    best=None
    for v in freevars:
        old=val[v]; val[v]=old+1; d=evn(node)-cur; val[v]=old
        if d!=0 and cur%d==0:
            # prefer free inputs not already deliberately set, and appearing in fewer eqs
            score=(v in setfree, abs(cur//d).bit_length())
            if best is None or score<best[0]: best=(score,v,d)
    if best is None: return False
    _,v,d=best
    val[v]-=cur//d
    return True
forward(); ns['v']=val
prev=None
for it in range(200):
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if it%10==0 or len(F)<30: print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
    if not F: print("SOLVED!"); break
    changed=False
    for i in F:
        for t in flat(rootast(i)):
            if evn(t)!=0:
                if try_absorb(t): changed=True; break
    forward()
    if not changed:
        print(f"  no absorber found; stuck at {len(F)} fail"); break
    if prev is not None and len(F)>=prev and it>5:
        pass  # allow non-monotone briefly
    prev=len(F)
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:20]}", flush=True)
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('constructive_solved.json','w')); print("SAVED")
elif len(F)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('constructive_partial.json','w')); print("saved partial")
