#!/usr/bin/env python3
"""Global constrained solve from the ALL-ZERO (0,0) point (x_15298=0, all free=0). Linearization is
exact (all products vanish). Target: the failing equation roots -> 0 while keeping satisfied roots
at 0. Handles = closure of the failing set. If consistent, apply + Newton; else report definitively."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
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
posof={t:k for k,t in enumerate(order)}
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
consumers=defaultdict(list)
for k,t in enumerate(order):
    for u in gates[definer[t]][2]: consumers[u].append(k)
def downstream_ks(w):
    seenv=set(); affected=set(); stack=[w]
    while stack:
        x=stack.pop()
        if x in seenv: continue
        seenv.add(x)
        for k in consumers.get(x,()):
            if k not in affected: affected.add(k); stack.append(order[k])
    return sorted(affected)
ns={'__builtins__':{}}; ns['v']=val
def forward():
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
forward()
byh=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]&freeinp: byh[v].add(i)
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"all-zero (0,0): {len(lines)-len(F0)}/{len(lines)} ({len(F0)} fail), x_15298={val[15298]}", flush=True)
# handle closure of F0 (2 hops)
W=set()
for i in F0: W|=(eqvars[i]&freeinp)
for _ in range(2):
    eqs=set()
    for h in list(W): eqs|=byh[h]
    for i in eqs: W|=(eqvars[i]&freeinp)
W=sorted(W-set(v for v in freeinp if pinned[v]))
print(f"closure handles: {len(W)}", flush=True)
dks={}; deq={}
for w in W:
    ks=downstream_ks(w); dks[w]=ks
    aff=set(eqbyvar.get(w,set()))
    for k in ks: aff|=eqbyvar.get(order[k],set())
    deq[w]=sorted(aff)
def inv(a): return pow(a%p,p-2,p)
def rref_solve(rows, rhs, ncol):
    M=[[rows[i][j]%p for j in range(ncol)]+[rhs[i]%p] for i in range(len(rows))]
    m=len(M); r=0; where=[-1]*ncol
    for c in range(ncol):
        piv=None
        for i in range(r,m):
            if M[i][c]%p: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; iv=inv(M[r][c]); M[r]=[(x*iv)%p for x in M[r]]
        for i in range(m):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][k]-f*M[r][k])%p for k in range(ncol+1)]
        where[c]=r; r+=1
        if r>=m: break
    for i in range(r,m):
        if M[i][ncol]%p!=0: return None,r
    x=[0]*ncol
    for c in range(ncol):
        if where[c]!=-1: x[c]=M[where[c]][ncol]%p
    return x,r
for newton in range(6):
    forward()
    baseroot={i:eval(rootcode_of(i),ns)%p for i in range(len(lines))}
    Fn=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if not Fn: print(f"newton {newton}: SOLVED (mod p)"); break
    cols=[]
    for w in W:
        o=val[w]; val[w]=o+1
        for k in dks[w]: val[order[k]]=eval(gcode[k],ns)
        col={}
        for i in deq[w]:
            d=(eval(rootcode_of(i),ns)-baseroot[i])%p
            if d: col[i]=d
        cols.append(col); val[w]=o
        for k in dks[w]: val[order[k]]=eval(gcode[k],ns)
    touched=set()
    for col in cols: touched|=set(col)
    rows=[]; rhs=[]
    for i in touched:
        rows.append([cols[j].get(i,0) for j in range(len(W))])
        rhs.append((-baseroot[i])%p)   # want root -> 0 (fixes fails, keeps satisfied at 0)
    x,rk=rref_solve(rows,rhs,len(W))
    print(f"newton {newton}: |eqs|={len(touched)}, |W|={len(W)}, rank={rk}, {'CONSISTENT' if x is not None else 'INCONSISTENT'}", flush=True)
    if x is None: break
    for j,w in enumerate(W):
        d=x[j]%p
        if d>p//2: d-=p
        val[w]+=d
    forward()
forward()
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('q00_GLOBAL_SOLVED.json','w')); print("*** SOLVED ***")
