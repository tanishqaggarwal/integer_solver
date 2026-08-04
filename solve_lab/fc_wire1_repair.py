#!/usr/bin/env python3
"""Local repair on the wire=1 branch (methodology escape). Load best, pin wire=sign*1, forward-eval,
then Dixon-repair the failing equations (13 unpackings + core, which is trivial with wire=1)."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
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
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v] and v not in wire: val[v]=(-c0)//c1; pinned[v]=True
WIREV=1
for v,s in wire.items(): val[v]=s*WIREV; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
for v,x in override.items():
    if v not in wire: val[v]=x; pinned[v]=True
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
    if v in best and v not in wire: val[v]=best[v]
forward(); ns['v']=val
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"best + wire=1: {len(lines)-len(F0)}/{len(lines)} ({len(F0)} fail)", flush=True)
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
def dixon(M,b,steps=14):
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
Faccum=set()
for rnd in range(20):
    ns['v']=val
    Fnow=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if not Fnow: print("SOLVED!"); break
    Faccum|=set(Fnow); F=sorted(Faccum)
    H=sorted(set().union(*[eqvars[i]&freeinp for i in F])-set(wire))
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
    print(f"rnd {rnd}: {len(Fnow)} now, {len(F)} accum, {len(H)} handles, rank {r}", flush=True)
    if y is None: print("  dixon failed"); break
    for j in range(r): val[H[pc[j]]]+=y[j]
    forward()
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {sorted(F)[:20]}")
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('wire1_SOLVED.json','w')); print("*** SAVED wire1_SOLVED.json ***")
