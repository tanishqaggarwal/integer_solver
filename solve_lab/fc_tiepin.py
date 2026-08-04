#!/usr/bin/env python3
"""Forward-construct with difference-square partners PINNED equal to their load, then re-solve the
linear equations (constructor re-absorbs). Kills the (load-partner)^2 monsters at the source."""
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
gate_defs=defaultdict(list)
for t,rhs,vids in gates: gate_defs[t].append(rhs)
# find difference-square pairs: gate y=z*z (z=z) with z = a - b; collect (a,b) where b is a free input
sol0={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
pairpins={}
for t in gate_out:
    rhs=gate_defs[t][0]
    m=re.match(r'x_(\d+) \* x_(\d+)$', rhs)
    if m and m.group(1)==m.group(2):
        z=int(m.group(1))
        zd=gate_defs.get(z,[''])[0]
        md=re.match(r'x_(\d+) - x_(\d+)$', zd)
        if md:
            a,b=int(md.group(1)),int(md.group(2))
            av=sol0.get(a,0)
            if abs(av)>2**200:  # a is a load
                if b in freeinp: pairpins[b]=av
                elif a in freeinp: pairpins[a]=sol0.get(b,0)
print(f"difference-square partner pins: {len(pairpins)}", flush=True)
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
override.update(pairpins)
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
anc=defaultdict(set)
for v in freeinp: anc[v]={v}
for k,t in enumerate(order):
    s=set()
    for u in gates[definer[t]][2]: s|=anc[u]
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
def dep_final(w): return w in determined or val[w]==0
def try_set(t):
    frees=[v for v in gvars(t) if v in freeinp and v not in determined]
    if not frees: return False
    rem=None; quot=None
    for v in frees:
        if not all(dep_final(w) for w in free_deps(t)-{v}): continue
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
            dr=(-cur*pow(cr%p,p-2,p))%p; Rp=cur+cr*dr
            if Rp%cq==0: val[vr]+=dr; val[vq]+=(-Rp)//cq; determined.add(vr); determined.add(vq); return True
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
    if it%40==0 or len(F)<30: print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
    if not F: print("SOLVED!"); break
    changed=False
    for i in F:
        for t in flat(rootast(i)):
            if evn(t)!=0 and try_set(t): changed=True
    forward()
    if not changed: print(f"  stuck at {len(F)} fail"); break
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:20]}", flush=True)
if len(F)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('fc_tiepin_solved.json','w')); print("SAVED")
