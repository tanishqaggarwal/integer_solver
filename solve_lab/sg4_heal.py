import heal_harness as H
import ast, re, sys
from collections import defaultdict
p=H.p
# ---- build equation flat-term machinery (like forward_construct) ----
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
val=H.val
ns={'v':val,'__builtins__':{}}
def evn(node):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -evn(node.operand)
    a=evn(node.left); b=evn(node.right)
    return a+b if isinstance(node.op,ast.Add) else a-b if isinstance(node.op,ast.Sub) else a*b
def flat(node,s=1,o=None):
    if o is None:o=[]
    if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub)):
        flat(node.left,s,o); flat(node.right,s*(1 if isinstance(node.op,ast.Add) else -1),o)
    else: o.append((node,s))
    return o
astcache={}
def rootterms(i):
    if i not in astcache:
        node=ast.parse(lines[i].rsplit('=',1)[0],mode='eval').body
        astcache[i]=flat(node)
    return astcache[i]
def gvars(node): return set(int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node)))
def coeff(node,v):
    old=val[v]; base=evn(node); val[v]=old+1; c=evn(node)-base; val[v]=old; return c
# free-input ancestor deps
def free_deps(node):
    s=set()
    for v in gvars(node): s|=H.anc.get(v,{v} if v in H.freeinp else set())
    return s

def load_x4287_state():
    vA=H.loadd('best_agentA_39022.json')
    for v in H.freeinp: val[v]=vA.get(v,0)
    val[4287]=1; val[9413]=0; val[17325]=0
    H.forward()
    # post-flip gaps
    x8731=val[4432]-val[20492]; x9118=val[7068]-val[37158]
    val[8731]=x8731; val[9118]=x9118
    H.forward()

load_x4287_state()
F=set(H.fails())
print(f"x_4287=1 start: {len(F)} fails")
PROTECT={4287,8731,9118,9413,17325,4432,7068,14853,12186,16742,2081,24601,30213,22162,24468,18956}
determined=set(PROTECT)
def dep_final(w): return w in determined or val[w]==0
def try_set(term_node,sgn):
    frees=[v for v in gvars(term_node) if v in H.freeinp and v not in determined]
    if not frees: return False
    rem=None; quot=None
    fd=free_deps(term_node)
    others_det=lambda v: all(dep_final(w) for w in fd-{v})
    cur_full=None
    for v in frees:
        if not others_det(v): continue
        c=coeff(term_node,v)  # coeff within the term
        if c==0: continue
        if c%p==0:
            if quot is None: quot=(v,c)
        else:
            if rem is None: rem=(v,c)
    return rem,quot

# heal loop: for each failing eq, try to zero it by setting a free input (double-width)
import time
t0=time.time()
best=len(F)
for it in range(300):
    H.forward()
    F=[i for i in range(len(lines)) if eval(H.eqcode[i],ns)!=0]
    if len(F)<best:
        best=len(F)
    if it%10==0: print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail) t={time.time()-t0:.0f}s",flush=True)
    if not F: print("SOLVED!"); break
    changed=False
    for i in F:
        # whole-equation residual
        node=ast.parse(lines[i].rsplit('=',1)[0],mode='eval').body
        cur=evn(node)
        if cur==0: continue
        # find a free input in this eq to zero it (double-width split over whole eq)
        eqfree=[v for v in gvars(node) if v in H.freeinp and v not in determined]
        rem=None;quot=None
        for v in eqfree:
            c=coeff(node,v)
            if c==0: continue
            if c%p==0:
                if quot is None: quot=(v,c)
            elif rem is None: rem=(v,c)
        if rem is not None and quot is not None:
            vr,cr=rem; vq,cq=quot
            dr=(-cur*pow(cr%p,p-2,p))%p
            Rp=cur+cr*dr
            if Rp%cq==0:
                val[vr]+=dr; val[vq]+=(-Rp)//cq; determined.add(vr);determined.add(vq);changed=True; continue
        if rem is not None:
            vr,cr=rem
            if cur%cr==0: val[vr]-=cur//cr; determined.add(vr); changed=True; continue
        if quot is not None:
            vq,cq=quot
            if cur%cq==0: val[vq]-=cur//cq; determined.add(vq); changed=True; continue
    if not changed:
        print(f"  stuck at {len(F)} fail: {F[:20]}"); break
H.forward()
F=[i for i in range(len(lines)) if eval(H.eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:25]}")
import json
if len(F)<11:
    json.dump({f"x_{i}":val[i] for i in range(H.NVARS)},open('sg4_heal_out.json','w'))
    print("saved sg4_heal_out.json")
