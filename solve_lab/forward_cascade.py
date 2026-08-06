#!/usr/bin/env python3
"""Cascade closer: start from agentA (core solved, S=T=0), PIN the core cone as determined,
then run forward_construct's double-width try_set heal loop to close the G1/G2 cascade around
the fixed core. Process in gate-topological order to avoid oscillation."""
import json, re, ast, sys
from collections import defaultdict, deque
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
import heal_harness as H
val=H.val; freeinp=H.freeinp; order=H.order; anc=H.anc; gcode=H.gcode
lines=H.lines; eqcode=H.eqcode; eqvars=H.eqvars; NVARS=H.NVARS
gates=H.gates; definer_t={t:i for i,(t,_,_) in enumerate(gates)}
# map target -> its defining gate (from heal_harness order/definer)
definer={}
gu=[0]*len(gates); ready=[False]*NVARS
for v in range(NVARS):
    if v in freeinp: ready[v]=True
using=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates):
    for v in vids: using[v].append(gi)
    gu[gi]=sum(1 for v in vids if not ready[v])
q=deque(gi for gi in range(len(gates)) if gu[gi]==0)
while q:
    gi=q.popleft(); t,rhs,vids=gates[gi]
    if ready[t]: continue
    definer[t]=gi; ready[t]=True
    for gj in using[t]:
        gu[gj]-=1
        if gu[gj]==0: q.append(gj)

# load agentA
vA=H.loadd('best_agentA_39022.json')
for v in freeinp: val[v]=vA.get(v,0)
H.forward()
ns={'v':val,'__builtins__':{}}
def forward(): H.forward()
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
def gvars(node): return set(int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node)))
def free_deps(node):
    s=set()
    for v in gvars(node): s|=anc.get(v,{v} if v in freeinp else set())
    return s
def coeff(node, v):
    old=val[v]; base=evn(node); val[v]=old+1; c=evn(node)-base; val[v]=old; return c

# PIN core cone
core_cone=set()
for g in [35389,6671,24908,29322,3558,29356,27762,33469]:
    core_cone|=anc.get(g,set())
core_cone|={14853,12186,16742,30317,2936,5146}
# also message/quadrant pins
override={24601,2081,30213,22162,24468,18956,4287}
determined=set(v for v in freeinp if v in core_cone or v in override)
print("pinned core-cone free inputs:",len(determined))

def dep_final(w): return w in determined or val[w]==0
def try_set(t):
    frees=[v for v in gvars(t) if v in freeinp and v not in determined]
    if not frees: return False
    rem=None; quot=None
    others_det=lambda v: all(dep_final(w) for w in free_deps(t)-{v})
    for v in frees:
        if not others_det(v): continue
        c=coeff(t,v)
        if c==0: continue
        if c%p==0:
            if quot is None: quot=(v,c)
        else:
            if rem is None: rem=(v,c)
    cur=evn(t)
    if cur==0: return False
    if rem is not None and quot is not None:
        vr,cr=rem; vq,cq=quot
        if all(dep_final(w) for w in free_deps(t)-{vr,vq}):
            dr=(-cur*pow(cr%p,p-2,p))%p
            Rp=cur+cr*dr
            if Rp%cq==0:
                val[vr]+=dr; val[vq]+=(-Rp)//cq; determined.add(vr); determined.add(vq); return True
    if rem is not None:
        vr,cr=rem
        if cur%cr==0: val[vr]-=cur//cr; determined.add(vr); return True
    if quot is not None:
        vq,cq=quot
        if cur%cq==0: val[vq]-=cur//cq; determined.add(vq); return True
    return False

forward()
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print("start fails:",len(F0), sorted(F0))
for it in range(200):
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if it%10==0 or len(F)<25: print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); det={len(determined)}", flush=True)
    if not F: print("SOLVED!"); break
    changed=False
    for i in F:
        for t in flat(rootast(i)):
            if evn(t)!=0 and try_set(t): changed=True
    forward()
    if not changed: print(f"  stuck at {len(F)} fail: {sorted(F)[:25]}"); break
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {sorted(F)[:25]}")
if len(F)<11:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('/home/user/integer_solver/solve_lab/cascade_out.json','w'))
    print("saved cascade_out.json")
