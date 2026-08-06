#!/usr/bin/env python3
"""Quadrant (0,0) full solve: greedy-construct, then Dixon-solve the residual failures (which reduce
to a small set of directly-settable gadgets since x_15298=0). Verify + save."""
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
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
override={}  # quadrant (0,0): x_15298=0
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
def inner_src(lhs):
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
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
_rc={}
def rootcode_of(i):
    if i not in _rc: _rc[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rc[i]
def evn(node):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -evn(node.operand)
    a=evn(node.left); b=evn(node.right)
    return a+b if isinstance(node.op,ast.Add) else a-b if isinstance(node.op,ast.Sub) else a*b
def inner(lhs): return inner_src(lhs)
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
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
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
# greedy phase
forward()
for it in range(200):
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
print(f"after greedy: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); x_15298={val[15298]}", flush=True)
# Dixon phase on residual
def inv(a): return pow(a%p,p-2,p)
def gfp_pivots(Jac):
    m=len(Jac); n=len(Jac[0]) if m else 0
    M=[[Jac[i][j]%p for j in range(n)] for i in range(m)]
    rowmap=list(range(m)); pr=0; pivrows=[]; pivcols=[]
    for c in range(n):
        piv=None
        for i in range(pr,m):
            if M[i][c]%p!=0: piv=i;break
        if piv is None: continue
        M[pr],M[piv]=M[piv],M[pr]; rowmap[pr],rowmap[piv]=rowmap[piv],rowmap[pr]
        iv=inv(M[pr][c]); M[pr]=[(x*iv)%p for x in M[pr]]
        for i in range(m):
            if i!=pr and M[i][c]%p!=0:
                f=M[i][c]; M[i]=[(M[i][k]-f*M[pr][k])%p for k in range(n)]
        pivrows.append(rowmap[pr]); pivcols.append(c); pr+=1
        if pr>=m: break
    return pivrows,pivcols
def matinv(M):
    r=len(M); Aug=[[M[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv=None
        for i in range(c,r):
            if Aug[i][c]%p!=0: piv=i;break
        if piv is None: return None
        Aug[c],Aug[piv]=Aug[piv],Aug[c]; iv=inv(Aug[c][c]); Aug[c]=[(x*iv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p!=0:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(M,b,steps=12):
    r=len(M); Mi=matinv(M)
    if Mi is None: return None
    x=[0]*r; bb=b[:]; mod=1
    for _ in range(steps):
        bm=[bb[i]%p for i in range(r)]
        xi=[sum(Mi[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        nb=[]
        for i in range(r):
            s=bb[i]-sum(M[i][k]*xi[k] for k in range(r))
            if s%p!=0: return None
            nb.append(s//p)
        bb=nb; mod*=p
        if all(z==0 for z in bb): break
    half=mod//2; y=[]
    for xi in x:
        xi%=mod
        if xi>half: xi-=mod
        y.append(xi)
    return y
for rnd in range(15):
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if not F: print("SOLVED via Dixon!"); break
    H=sorted(set().union(*[eqvars[i]&freeinp for i in F]))
    base=[eval(rootcode_of(i),ns) for i in F]
    Jac=[[0]*len(H) for _ in F]
    for j,h in enumerate(H):
        o=val[h]; val[h]=o+1; forward(); ns['v']=val
        for ri in range(len(F)): Jac[ri][j]=eval(rootcode_of(F[ri]),ns)-base[ri]
        val[h]=o
    forward(); ns['v']=val
    pr,pc=gfp_pivots(Jac); r=len(pr)
    M=[[Jac[pr[i]][pc[j]] for j in range(r)] for i in range(r)]
    rhs=[-base[pr[i]] for i in range(r)]
    y=dixon(M,rhs)
    print(f"rnd {rnd}: {len(F)} fail, {len(H)} handles, rank {r}", flush=True)
    if y is None: print("  dixon failed"); break
    for j in range(r): val[H[pc[j]]]+=y[j]
    forward()
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:20]}; x_15298={val[15298]}")
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('q00_SOLVED.json','w')); print("*** SAVED q00_SOLVED.json ***")
