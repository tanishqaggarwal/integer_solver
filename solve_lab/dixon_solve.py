#!/usr/bin/env python3
"""Iterative integer solve using Dixon p-adic lifting (fast for huge-entry systems).
Quadrant (1,1) config. Each round: accumulate failing eqs, build Jacobian, select full-rank
square subsystem via GF(p) elimination, Dixon-solve, apply, forward-eval. Verify."""
import json, re, ast, sys
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
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
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
_rootcache={}
def rootcode_of(i):
    if i not in _rootcache:
        _rootcache[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rootcache[i]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def inv_modp(a): return pow(a%p, p-2, p)
def matinv_modp(M):
    r=len(M); Aug=[[M[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv=None
        for i in range(c,r):
            if Aug[i][c]%p!=0: piv=i;break
        if piv is None: return None
        Aug[c],Aug[piv]=Aug[piv],Aug[c]
        iv=inv_modp(Aug[c][c]); Aug[c]=[(x*iv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p!=0:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(M, b, steps=6):
    """Solve M y = b over Z (M nonsingular mod p). Returns y or None."""
    r=len(M); Minv=matinv_modp(M)
    if Minv is None: return None
    x=[0]*r; bb=b[:]; mod=1
    for _ in range(steps):
        # xi = Minv (bb mod p)
        bm=[bb[i]%p for i in range(r)]
        xi=[sum(Minv[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        # bb = (bb - M xi)/p
        nb=[]
        for i in range(r):
            s=bb[i]-sum(M[i][k]*xi[k] for k in range(r))
            if s%p!=0: return None
            nb.append(s//p)
        bb=nb; mod*=p
        if all(z==0 for z in bb): break
    # balanced representation mod 'mod'
    half=mod//2
    y=[]
    for xi in x:
        xi%=mod
        if xi>half: xi-=mod
        y.append(xi)
    return y
def gfp_pivots(Jac):
    """GF(p) elimination to find pivot rows/cols of Jac (m x n)."""
    m=len(Jac); n=len(Jac[0]) if m else 0
    M=[[Jac[i][j]%p for j in range(n)] for i in range(m)]
    rowmap=list(range(m)); pr=0; pivrows=[]; pivcols=[]
    for c in range(n):
        piv=None
        for i in range(pr,m):
            if M[i][c]%p!=0: piv=i;break
        if piv is None: continue
        M[pr],M[piv]=M[piv],M[pr]; rowmap[pr],rowmap[piv]=rowmap[piv],rowmap[pr]
        iv=inv_modp(M[pr][c])
        M[pr]=[(x*iv)%p for x in M[pr]]
        for i in range(m):
            if i!=pr and M[i][c]%p!=0:
                f=M[i][c]; M[i]=[(M[i][k]-f*M[pr][k])%p for k in range(n)]
        pivrows.append(rowmap[pr]); pivcols.append(c); pr+=1
        if pr>=m: break
    return pivrows, pivcols
def solve_step(H, F):
    n=len(H); m=len(F)
    ns['v']=val
    base=[eval(rootcode_of(i),ns) for i in F]
    Jac=[[0]*n for _ in range(m)]
    for j,h in enumerate(H):
        old=val[h]; val[h]=old+1; forward(); ns['v']=val
        for ri in range(m): Jac[ri][j]=eval(rootcode_of(F[ri]),ns)-base[ri]
        val[h]=old
    forward()
    pivrows,pivcols=gfp_pivots(Jac)
    r=len(pivrows)
    M=[[Jac[pivrows[i]][pivcols[j]] for j in range(r)] for i in range(r)]
    rhs=[-base[pivrows[i]] for i in range(r)]
    y=dixon(M,rhs)
    if y is None: return False,r
    for j in range(r): val[H[pivcols[j]]]+=y[j]
    forward()
    return True,r
forward(); ns['v']=val
Faccum=set()
for it in range(25):
    ns['v']=val
    Ffail=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    print(f"iter {it}: {len(lines)-len(Ffail)}/{len(lines)} ({len(Ffail)} fail)", flush=True)
    if not Ffail: print("SOLVED!"); break
    Faccum |= set(Ffail); Flist=sorted(Faccum)
    H=sorted(set().union(*[eqvars[i]&hbase for i in Flist]))
    print(f"  system: {len(Flist)} eqs, {len(H)} handles", flush=True)
    ok,r=solve_step(H, Flist)
    if not ok: print(f"  Dixon failed (rank {r})"); break
ns['v']=val
Ffail=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(Ffail)}/{len(lines)} ({len(Ffail)} fail): {Ffail[:20]}", flush=True)
if len(Ffail)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('dixon_solved.json','w')); print("SAVED dixon_solved.json")
elif len(Ffail)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('dixon_partial.json','w')); print("saved partial")
