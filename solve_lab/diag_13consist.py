#!/usr/bin/env python3
"""From best + wire=1, check if the 13 breaking unpackings are CONSISTENT (solvable) via their
free-input handles (excl wire). Jacobian rank vs augmented rank. Consistent => wire escape closes."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
A=load_atoms()
par=list(range(NVARS)); sgn=[1]*NVARS
def find2(x):
    s=1; r=x
    while par[r]!=r: s*=sgn[r]; r=par[r]
    return r,s
def union(a,b,rel):
    ra,sa=find2(a); rb,sb=find2(b)
    if ra==rb: return
    par[ra]=rb; sgn[ra]=rel*sb*sa
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==2 and pp.get((),0)==0:
        v1,v2=sorted(vs); c1=pp.get((v1,),0); c2=pp.get((v2,),0)
        qok=all(pp.get(k,0)==0 for k in pp if isinstance(k,tuple) and len(k)==2)
        if qok and c1!=0 and c2!=0 and abs(c1)==abs(c2):
            rel=(-c2)//c1
            if rel in (1,-1): union(v1,v2,rel)
r0,_=find2(26064); wire={v:find2(v)[1] for v in range(NVARS) if find2(v)[0]==r0}
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v] and v not in wire: val[v]=(-c0)//c1; pinned[v]=True
for v,s in wire.items(): val[v]=s*1; pinned[v]=True  # wire=1
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
    if ready[t] or t in wire: continue
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
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
_rc={}
def rc(i):
    if i not in _rc: _rc[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rc[i]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best and v not in wire: val[v]=best[v]
forward(); ns['v']=val
b13=[8429, 11166, 11915, 12594, 23869, 25313, 26785, 31400, 32300, 36106, 36767, 37257]
b13=[i for i in b13 if eval(compile(VAR.sub(r'v[\1]',lines[i].rsplit('=',1)[0]),'<e>','eval'),ns)!=0]
H=sorted(set().union(*[eqvars[i]&freeinp for i in b13])-set(wire))
print(f"wire=1: {len(b13)} of the 13 broken; {len(H)} handles")
base=[eval(rc(i),ns) for i in b13]
Jac=[[0]*len(H) for _ in b13]
for j,h in enumerate(H):
    o=val[h]; val[h]=o+1; forward(); ns['v']=val
    for ri in range(len(b13)): Jac[ri][j]=eval(rc(b13[ri]),ns)-base[ri]
    val[h]=o
forward()
def rank_modp(M):
    M=[[x%p for x in row] for row in M]; m=len(M); n=len(M[0]) if m else 0; r=0
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
aug=[Jac[i]+[-base[i]] for i in range(len(b13))]
rka=rank_modp(aug)
print(f"13-unpacking Jacobian rank={rk}/{len(b13)}, augmented={rka} => {'CONSISTENT (solvable!)' if rka==rk else 'INCONSISTENT'}")
