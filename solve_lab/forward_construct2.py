#!/usr/bin/env python3
"""Topological forward-construction. Precompute each var's free-input ancestors. Set each free
input to its gadget value ONLY when all its dependencies are determined (value is final -> no
oscillation). Double-width split (r=target mod p, q=target//p) at wire-product chains. Verify."""
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
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
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
# precompute free-input ancestors of each gate target (its backward-cone leaves), in topo order
anc=defaultdict(set)  # var -> set of free inputs it depends on
for v in freeinp: anc[v]={v}
for k,t in enumerate(order):
    _,rhs,vids=gates[definer[t]]
    s=set()
    for u in vids: s|=anc[u]
    anc[t]=s
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
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
def gvars(node): return set(int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node)))
def free_deps(node):
    s=set()
    for v in gvars(node): s|=anc.get(v,{v} if v in freeinp else set())
    return s
def coeff(node, v):
    old=val[v]; base=evn(node); val[v]=old+1; c=evn(node)-base; val[v]=old; return c
determined=set(v for v in freeinp if pinned[v])
# ===== BAKE IN CORE (follower style) =====
# core: x_29322 = x_14853 - x_12186 = 0 ; x_3558 = x_24908 - x_16742 = 0
# x_12186 FOLLOWS x_14853 (so x_29322=0 exactly); x_16742 FOLLOWS x_24908 (so x_3558=0 exactly).
# The followers block (dep_final False) until their target is final, so the sweep uses OTHER
# rem-handles for the gadgets that x_12186/x_16742 vacate.
for cv in (30317,2936,5146): val[cv]=0; determined.add(cv)  # L=0 -> core quot handles 0
PENDING={12186,16742}
FOLLOW={12186:14853, 16742:24908}
FAILVARS=set()  # free inputs in currently-failing equations (may still change)
# a gadget-term t can determine free input v iff: v in freeinp, not determined, and the term's OTHER
# free-input dependencies are all FINAL (determined, or 0-and-not-in-a-failing-eq). Then set v.
def dep_final(w):
    if w in PENDING: return False
    return w in determined or val[w]==0
def target_final(t):
    # free target: final iff determined; gate target: all free ancestors determined-or-zero
    if t in freeinp: return t in determined
    return all((a in determined) or (val[a]==0) for a in anc.get(t,set()))
def resolve_followers():
    ch=False
    for f,t in list(FOLLOW.items()):
        if f not in PENDING: continue
        if target_final(t):
            forward(); val[f]=val[t]; determined.add(f); PENDING.discard(f); ch=True
    return ch
def try_set(t):
    frees=[v for v in gvars(t) if v in freeinp and v not in determined and v not in PENDING]
    if not frees: return False
    # candidate: a free input v s.t. all OTHER free-input deps of t are FINAL, and v has nonzero coeff
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
    # double-width split: rem (coeff invertible mod p) + quot (coeff mult of p), both deps-ready
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
forward(); ns['v']=val
for it in range(400):
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    FAILVARS=set()
    for i in F: FAILVARS|=(eqvars[i]&freeinp)
    FAILVARS-=determined
    if it%20==0 or len(F)<40: print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); determined free={len(determined)}", flush=True)
    if not F: print("SOLVED!"); break
    changed=False
    for i in F:
        for t in flat(rootast(i)):
            if evn(t)!=0 and try_set(t): changed=True
    if not changed:
        # sweep stuck: resolve deferred core followers (x_12186<-x_14853, x_16742<-x_24908)
        # once their targets are final, then let the sweep continue to heal downstream.
        if PENDING:
            for f,t in list(FOLLOW.items()):
                if f in PENDING: forward(); val[f]=val[t]; determined.add(f); PENDING.discard(f)
            forward(); changed=True; continue
        print(f"  no ready free input to set; stuck at {len(F)} fail"); break
    forward()
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:20]}", flush=True)
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('fc_solved.json','w')); print("SAVED fc_solved.json")
elif len(F)<30:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('fc_partial.json','w')); print("saved partial")
