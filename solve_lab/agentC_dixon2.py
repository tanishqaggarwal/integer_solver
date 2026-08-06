#!/usr/bin/env python3
"""Single-shot Dixon heal: from agentA_39021 + loads=CONST (24 broken), solve for non-load-cone
handle deltas that zero the 24 broken eqs' inner-roots while preserving all satisfied eqs touched.
No accumulation: build Jacobian once, Dixon-solve the rank-r pivot system once, apply, verify."""
import json, re, ast, sys, time
from collections import defaultdict, deque
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone, load_solution,
                            forward, eqcode, eqvars, lines, NEQ, NVARS, anc)
sys.setrecursionlimit(1000000)
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
CORE={2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
QUOT={30317,2936,5146}
gvids={t:gates[definer[t]][2] for t in order}
posof={t:k for k,t in enumerate(order)}
ns={'__builtins__':{}}
VAR=re.compile(r'x_(\d+)')
def count(v):
    ns['v']=v; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)
def set_quot(v):
    if v[11150]%p==0: v[30317]=-(v[11150])//p
    if (537773*v[37758])%p==0: v[2936]=(537773*v[37758])//p
    if v[25739]%(6672769*p)==0: v[5146]=v[25739]//(6672769*p)
def free_cone(r):
    seen=set(); st=[r]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u,()):
            if w not in seen: st.append(w)
    return set(u for u in seen if u in freeinp)
# inner root (strip outer const-mult and squares) like forward_construct
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
_rc={}
def rootcode(i):
    if i not in _rc:
        _rc[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rc[i]
# consumers for downstream
consumers=defaultdict(list)
for k,t in enumerate(order):
    for u in gvids[t]: consumers[u].append(k)
def downstream_ks(w):
    aff=set(); st=[w]; seen=set()
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        for k in consumers.get(x,()):
            if k not in aff: aff.add(k); st.append(order[k])
    return sorted(aff)
def partial_forward(v, ks):
    ns['v']=v
    for k in ks: v[order[k]]=eval(gcode_k[k],ns)
gcode_k=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
inv=lambda a: pow(a%p,p-2,p)

d=json.load(open('agentC_healdix.json'))
BROKEN=d['broken']; H=d['handles']
v=load_solution('best_agentA_39021.json'); v[22152]=CONST2; v[33462]=CONST1; forward(v); set_quot(v)
ns['v']=v
cur=count(v)
print(f"start: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail); broken={len(BROKEN)}, handles={len(H)}", flush=True)
eqbyvar=defaultdict(set)
for i in range(NEQ):
    for var in eqvars[i]: eqbyvar[var].add(i)
# constraint set: BROKEN (target 0) + satisfied eqs touched by any handle downstream (preserve 0)
dks={h:downstream_ks(h) for h in H}
touch=set(BROKEN)
for h in H:
    aff=set(eqbyvar.get(h,()))
    for k in dks[h]: aff|=eqbyvar.get(order[k],set())
    touch|=aff
sat_touch=sorted((touch-set(BROKEN)) - cur)  # currently-satisfied eqs touched by handles
cons_eqs=sorted(BROKEN)+sat_touch
print(f"constraint eqs: {len(BROKEN)} broken + {len(sat_touch)} satisfied-touched = {len(cons_eqs)}", flush=True)
base_root={i:eval(rootcode(i),ns)%p for i in cons_eqs}
# Jacobian mod p
t0=time.time()
Jcols=[]
consset=set(cons_eqs)
for hi,h in enumerate(H):
    o=v[h]; v[h]=o+1; partial_forward(v,dks[h]); ns['v']=v
    col={}
    aff=set(eqbyvar.get(h,()))&consset
    for k in dks[h]: aff|=(eqbyvar.get(order[k],set())&consset)
    for i in aff:
        dd=(eval(rootcode(i),ns)-base_root[i])%p
        if dd: col[i]=dd
    Jcols.append(col); v[h]=o; partial_forward(v,dks[h]); ns['v']=v
print(f"Jacobian built {time.time()-t0:.0f}s", flush=True)
# rows
rows=sorted(set().union(*[set(c) for c in Jcols])) if Jcols else []
ridx={r:i for i,r in enumerate(rows)}
M=len(rows)
Jm=[[0]*len(H) for _ in range(M)]
for hj,col in enumerate(Jcols):
    for r,val_ in col.items(): Jm[ridx[r]][hj]=val_%p
# target b: for BROKEN want root->0 so delta = -base_root; for sat want 0
bvec=[0]*M
for i in cons_eqs:
    if i in ridx:
        bvec[ridx[i]]=(-base_root[i])%p if i in set(BROKEN) else 0
# gaussian pivot to find rank & solve (mod p) then Dixon lift
def gfp_pivots(Jac):
    m=len(Jac); n=len(Jac[0]) if m else 0
    Mx=[[Jac[i][j]%p for j in range(n)] for i in range(m)]
    rowmap=list(range(m)); pr=0; pivrows=[]; pivcols=[]
    for c in range(n):
        piv=None
        for i in range(pr,m):
            if Mx[i][c]%p: piv=i;break
        if piv is None: continue
        Mx[pr],Mx[piv]=Mx[piv],Mx[pr]; rowmap[pr],rowmap[piv]=rowmap[piv],rowmap[pr]
        ivv=inv(Mx[pr][c]); Mx[pr]=[(x*ivv)%p for x in Mx[pr]]
        for i in range(m):
            if i!=pr and Mx[i][c]%p:
                f=Mx[i][c]; Mx[i]=[(Mx[i][k]-f*Mx[pr][k])%p for k in range(n)]
        pivrows.append(rowmap[pr]); pivcols.append(c); pr+=1
        if pr>=m: break
    return pivrows,pivcols
pr,pc=gfp_pivots(Jm)
rank=len(pr)
print(f"Jacobian rank={rank} over {len(H)} handles, {M} constraint-rows", flush=True)
# consistency: does bvec lie in column space? Solve pivot system, then verify all rows.
# Build square pivot system
Msq=[[Jm[pr[i]][pc[j]] for j in range(rank)] for i in range(rank)]
rhs=[bvec[pr[i]] for i in range(rank)]
def matinv(Mmat):
    r=len(Mmat); Aug=[[Mmat[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv=None
        for i in range(c,r):
            if Aug[i][c]%p: piv=i;break
        if piv is None: return None
        Aug[c],Aug[piv]=Aug[piv],Aug[c]; ivv=inv(Aug[c][c]); Aug[c]=[(x*ivv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(Mmat,b,steps=12):
    r=len(Mmat); Mi=matinv(Mmat)
    if Mi is None: return None
    x=[0]*r; bb=b[:]; mod=1
    for _ in range(steps):
        bm=[bb[i]%p for i in range(r)]
        xi=[sum(Mi[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        nb=[]
        for i in range(r):
            s=bb[i]-sum(Mmat[i][k]*xi[k] for k in range(r))
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
y=dixon(Msq,rhs)
if y is None:
    print("Dixon failed (pivot system not integer-solvable)"); sys.exit(0)
# verify pivot solution satisfies ALL constraint rows mod p
delta={H[pc[j]]:y[j] for j in range(rank)}
ok=True
for r in range(M):
    s=sum(Jm[r][pc[j]]*y[j] for j in range(rank))%p
    if s!=bvec[r]%p: ok=False; break
print(f"pivot solution satisfies all rows mod p: {ok}")
# apply integer deltas
for h,dd in delta.items(): v[h]=v[h]+dd
forward(v); set_quot(v); ns['v']=v
F=count(v)
core=F&CORE; nc=F-CORE
print(f"AFTER Dixon apply: {NEQ-len(F)}/{NEQ} ({len(F)} fail) core={len(core)} noncore={len(nc)}")
print(f"  noncore: {sorted(nc)}")
if len(F)==0:
    json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open('best_agentC_39033.json','w'))
    print("*** FULL WIN best_agentC_39033.json ***")
elif NEQ-len(F)>39021:
    json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open(f'best_agentC_{NEQ-len(F)}.json','w'))
    print(f"IMPROVED -> best_agentC_{NEQ-len(F)}.json")
