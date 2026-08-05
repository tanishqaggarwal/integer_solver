#!/usr/bin/env python3
"""Reusable forward-construction harness for scanning configurations.
Preloads atoms/gates/eqcode/topo-order ONCE (override only touches free inputs, so the
gate topological order is config-independent). run_config(override) resets val/determined
and runs the greedy forward construction, returning metrics."""
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
val0=[0]*NVARS; pin0=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pin0[v]: val0[v]=(-c0)//c1; pin0[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand)
# --- topo order built ONCE using pin0 only (override adds only free inputs -> already ready) ---
ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp or pin0[v]: ready[v]=True
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
    _,rhs,vids=gates[definer[t]]
    s=set()
    for u in vids: s|=anc[u]
    anc[t]=s
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
CORE=[2071, 4573, 7123, 7469, 11854, 13660, 15299, 16622, 17726, 21382, 22093, 25480, 25539, 28653, 29437, 31061, 32894, 32916, 34517, 34892]
CORESET=set(CORE)
ns={'__builtins__':{}}
# global mutable state used by helpers
val=val0[:]
determined=set()
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
def dep_final(w):
    return w in determined or val[w]==0
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

def run_config(override, maxit=400, want_val=False):
    """Run forward construction with given override dict {var:value}. Returns metrics dict."""
    global val, determined
    val=val0[:]; pinned=pin0[:]
    for v,x in override.items(): val[v]=x; pinned[v]=True
    determined=set(v for v in freeinp if pinned[v])
    forward(); ns['v']=val
    for it in range(maxit):
        ns['v']=val
        F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
        if not F: break
        changed=False
        for i in F:
            for t in flat(rootast(i)):
                if evn(t)!=0 and try_set(t): changed=True
        forward()
        if not changed: break
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    Fset=set(F)
    core_fail=len(Fset & CORESET)
    noncore_fail=len(Fset - CORESET)
    S=val[35389]; T=val[6671]
    res={
        'satisfied': len(lines)-len(F),
        'nfail': len(F),
        'x_15298': val[15298],
        'x_7715': val[7715],
        'x_34554': val[34554],
        'S_modp': S%p, 'T_modp': T%p,
        'S_is0': (S%p==0), 'T_is0': (T%p==0),
        'core_fail': core_fail,
        'noncore_fail': noncore_fail,
        'F': F,
        'determined': len(determined),
    }
    if want_val: res['val']=val[:]
    return res

if __name__=='__main__':
    # sanity: reproduce baseline (1,1)
    r=run_config({24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2})
    print("baseline (1,1):", {k:v for k,v in r.items() if k!='F'})
    print("F:", r['F'])
