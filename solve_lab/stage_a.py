#!/usr/bin/env python3
"""STAGE A: solve the field (mod p) solution. All arithmetic mod p, so wire(=p) products vanish
and values stay < p (no huge-value ripple). Iteratively solve failing check-roots == 0 mod p for
the free-input remainders via GF(p) linear solve. Converges (bounded). Saves stage_a mod-p field."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1'])%p; C2=int(hc['C2'])%p
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
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=((-c0)//c1)%p; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
boolbits=set(json.load(open('boolbits.json'))['boolvars'])
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x%p; pinned[v]=True
hbase=freeinp - boolbits - set(override)
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
# gate code returns value mod p (wrap eval); we mod after each gate
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
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
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
_rc={}
def rc(i):
    if i not in _rc: _rc[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rc[i]
_ec={}
def ec(i):
    if i not in _ec: _ec[i]=compile(VAR.sub(r'v[\1]',lines[i].rsplit('=',1)[0]),'<e>','eval')
    return _ec[i]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)%p
def invp(a): return pow(a%p,p-2,p)
def gfp_solve(Jac, b):
    """Solve Jac y = b mod p (least, free=0). Jac m x n. Returns y (len n) or None."""
    m=len(Jac); n=len(Jac[0]) if m else 0
    Aug=[[Jac[i][j]%p for j in range(n)]+[b[i]%p] for i in range(m)]
    pr=0; piv=[]
    for c in range(n):
        sel=None
        for i in range(pr,m):
            if Aug[i][c]%p!=0: sel=i;break
        if sel is None: continue
        Aug[pr],Aug[sel]=Aug[sel],Aug[pr]
        iv=invp(Aug[pr][c]); Aug[pr]=[(x*iv)%p for x in Aug[pr]]
        for i in range(m):
            if i!=pr and Aug[i][c]%p!=0:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[pr][k])%p for k in range(n+1)]
        piv.append((pr,c)); pr+=1
        if pr>=m: break
    for i in range(m):
        if all(Aug[i][c]%p==0 for c in range(n)) and Aug[i][n]%p!=0: return None
    y=[0]*n
    for (r,c) in piv: y[c]=Aug[r][n]%p
    return y
def solve_step(H, F):
    n=len(H); m=len(F)
    ns['v']=val
    base=[eval(rc(i),ns)%p for i in F]
    Jac=[[0]*n for _ in range(m)]
    for j,h in enumerate(H):
        old=val[h]; val[h]=(old+1)%p; forward(); ns['v']=val
        for ri in range(m): Jac[ri][j]=(eval(rc(F[ri]),ns)-base[ri])%p
        val[h]=old
    forward()
    y=gfp_solve(Jac,[(-b)%p for b in base])
    if y is None: return False
    for j in range(n): val[H[j]]=(val[H[j]]+y[j])%p
    forward()
    return True
forward(); ns['v']=val
Faccum=set()
for it in range(40):
    ns['v']=val
    Ffail=[i for i in range(len(lines)) if eval(ec(i),ns)%p!=0]
    print(f"iter {it}: {len(lines)-len(Ffail)}/{len(lines)} fail-mod-p={len(Ffail)}", flush=True)
    if not Ffail: print("FIELD SOLUTION (all checks == 0 mod p)!"); break
    Faccum|=set(Ffail); Flist=sorted(Faccum)
    H=sorted(set().union(*[eqvars[i]&hbase for i in Flist]))
    print(f"  system {len(Flist)} eqs, {len(H)} handles", flush=True)
    if not solve_step(H, Flist):
        print("  GF(p) inconsistent (unexpected)"); break
ns['v']=val
Ffail=[i for i in range(len(lines)) if eval(ec(i),ns)%p!=0]
print(f"FINAL mod-p: {len(lines)-len(Ffail)}/{len(lines)} fail={len(Ffail)}", flush=True)
if len(Ffail)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('field_solution.json','w')); print("SAVED field_solution.json (mod-p field solution)")
