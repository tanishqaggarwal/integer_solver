#!/usr/bin/env python3
"""Decisive test: rank of the 23 broken roots' Jacobian wrt a LARGE handle set (all free inputs in
their 2-hop closure). Full rank(23) => healing possible with right handles; deficient => structural."""
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
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward()
val[14853]=val[14853]-(val[29322]%p); forward()
val[16742]=val[16742]+(val[3558]%p); forward(); ns['v']=val
broken=[3408, 3841, 4134, 4526, 5069, 7276, 15440, 15724, 15927, 21600, 22139, 22825, 27289, 27999, 28718, 29305, 31134, 31269, 32463, 33195, 36387, 36390, 38888]
byh=defaultdict(set)
for i,vs in enumerate(eqvars):
    for v in vs&freeinp: byh[v].add(i)
# 2-hop closure handle set
H=set()
for i in broken: H|=(eqvars[i]&freeinp)
for _ in range(2):
    eqs2=set()
    for h in list(H): eqs2|=byh[h]
    for i in eqs2: H|=(eqvars[i]&freeinp)
H=sorted(H-{14853,12186,16742,30317,2936,5146})
print(f"handles: {len(H)} for {len(broken)} broken")
base=[eval(rootcode_of(i),ns) for i in broken]
Jac=[[0]*len(H) for _ in range(len(broken))]
for j,h in enumerate(H):
    o=val[h]; val[h]=o+1; forward(); ns['v']=val
    for ri in range(len(broken)): Jac[ri][j]=eval(rootcode_of(broken[ri]),ns)-base[ri]
    val[h]=o
forward()
# rank mod p
def rank_modp(M):
    M=[[x%p for x in row] for row in M]; m=len(M); n=len(M[0]); r=0
    for c in range(n):
        piv=None
        for i in range(r,m):
            if M[i][c]%p: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; iv=pow(M[r][c],p-2,p); M[r]=[(x*iv)%p for x in M[r]]
        for i in range(m):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][k]-f*M[r][k])%p for k in range(n)]
        r+=1
        if r>=m: break
    return r
rk=rank_modp([row[:] for row in Jac])
print(f"rank of 23 broken wrt {len(H)} handles = {rk} / {len(broken)}")
# also check consistency: is base in column space? augment
aug=[Jac[i]+[-base[i]] for i in range(len(broken))]
rk_aug=rank_modp(aug)
print(f"augmented rank = {rk_aug}  => {'CONSISTENT (solvable)' if rk_aug==rk else 'INCONSISTENT'}")
