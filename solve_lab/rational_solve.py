#!/usr/bin/env python3
"""Solve the 27 build_twist failing equations EXACTLY over Q (Fraction) to distinguish true
inconsistency from integer p-divisibility. Report which handle gets a non-integer value."""
import json, re, ast, sys
from fractions import Fraction
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
boolbits=set(json.load(open('boolbits.json'))['boolvars'])
override={22106:1, 16742:C2, 12186:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
handles=sorted(freeinp - boolbits - set(override)); hset=set(handles)
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
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
forward(); ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"failing: {len(F)}", flush=True)
rootcode={i:compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval') for i in F}
H=sorted(set().union(*[eqvars[i]&hset for i in F]))
print(f"handles in failing eqs: {len(H)}", flush=True)
base={i:eval(rootcode[i],ns) for i in F}
# Jacobian (exact; roots linear in handles)
Jac=[[0]*len(H) for _ in F]
for j,h in enumerate(H):
    old=val[h]; val[h]=old+1; forward(); ns['v']=val
    for ri,i in enumerate(F): Jac[ri][j]=eval(rootcode[i],ns)-base[list(F).index(i)] if False else eval(rootcode[i],ns)-base[i]
    val[h]=old
forward()
# rational Gaussian elimination on Jac * x = -base
rows=[[Fraction(Jac[r][c]) for c in range(len(H))]+[Fraction(-base[F[r]])] for r in range(len(F))]
nr=len(rows); nc=len(H)
pr=0; piv=[]
for c in range(nc):
    sel=None
    for r in range(pr,nr):
        if rows[r][c]!=0: sel=r;break
    if sel is None: continue
    rows[pr],rows[sel]=rows[sel],rows[pr]
    pivval=rows[pr][c]; rows[pr]=[x/pivval for x in rows[pr]]
    for r in range(nr):
        if r!=pr and rows[r][c]!=0:
            f=rows[r][c]; rows[r]=[rows[r][k]-f*rows[pr][k] for k in range(nc+1)]
    piv.append(c); pr+=1
    if pr>=nr: break
# consistency over Q
incons=any(all(rows[r][k]==0 for k in range(nc)) and rows[r][nc]!=0 for r in range(nr))
print(f"Q-consistent: {not incons}  (rank {len(piv)} / {nc} handles, {len(F)} eqs)", flush=True)
if incons:
    print("TRULY INCONSISTENT over Q -> this activation/routing cannot work; try different config")
    sys.exit(0)
# particular solution (free vars=0): pivot vars = rows[r][nc]
sol={}
for r,c in enumerate(piv): sol[H[c]]=rows[r][nc]
noninteger=[(v,val_) for v,val_ in sol.items() if val_.denominator!=1]
print(f"Q-solution: {len(sol)} pivot handles set; NON-INTEGER handles: {len(noninteger)}", flush=True)
for v,fr in noninteger[:12]:
    print(f"  x_{v} = {fr.numerator}/{fr.denominator}  (denom factors: p? {fr.denominator%p==0}, denom={fr.denominator if fr.denominator<10**8 else '2^%d'%fr.denominator.bit_length()})")
json.dump({'noninteger':[(v,str(fr)) for v,fr in noninteger],'nfree':nc-len(piv)}, open('rational_sol.json','w'))
