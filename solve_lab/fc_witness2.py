#!/usr/bin/env python3
"""Sparse-witness solve. Only 30 slack free inputs are nonzero. Compute the Jacobian of all
currently-satisfied constraining equations + (x_3558, x_29322 residues) wrt these 30, mod p.
Find delta-w in the wiring null space that zeroes both residues. Apply (Newton), set quotients."""
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
forward(); ns['v']=val
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
# the 30 slack witness inputs (nonzero, non-override)
W0=set(v for v in freeinp if best.get(v,0)!=0 and v not in override)
# expand: free inputs in equations containing multi-role controls + x_24908 cone (movers)
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
def freecone(root):
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in cand:  # gate
            for u in gates[definer.get(x, [g for g in cand[x]][0] if False else cand[x][0])][2] if x in definer else []: st.append(u)
        elif x in freeinp: leaves.add(x)
    return leaves
extra=set()
for i in range(len(lines)):
    if eqvars[i]&{14853,12186,16742}:
        extra|=(eqvars[i]&freeinp)
# add x_24908 cone free inputs (movers of x_3558)
gd={t:vids for t,rhs,vids in gates}
def fcone(root):
    seen=set(); lv=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gd:
            for u in gd[x]: st.append(u)
        elif x in freeinp: lv.add(x)
    return lv
extra|=fcone(24908)
W=sorted((W0|extra)-set(override))
print(f"expanded witness inputs: {len(W)} (base {len(W0)})", flush=True)
def inv(a): return pow(a%p,p-2,p)
def rref_solve(rows, rhs, ncol):
    """Solve A x = rhs mod p (rows: list of ncol-vectors). Return (x or None, rank)."""
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
    # consistency
    for i in range(r,m):
        if M[i][ncol]%p!=0: return None,r
    x=[0]*ncol
    for c in range(ncol):
        if where[c]!=-1: x[c]=M[where[c]][ncol]%p
    return x,r
for newton in range(12):
    ns['v']=val
    # base roots for all equations
    baseroot={}
    for i in range(len(lines)):
        baseroot[i]=eval(rootcode_of(i),ns)
    b3=val[3558]%p; b29=val[29322]%p
    if b3==0 and b29==0: print(f"newton {newton}: residues ZERO"); break
    # Jacobian columns wrt W (mod p)
    cols=[]
    dres3=[]; dres29=[]
    for w in W:
        o=val[w]; val[w]=o+1; forward(); ns['v']=val
        col={}
        for i in range(len(lines)):
            d=(eval(rootcode_of(i),ns)-baseroot[i])%p
            if d: col[i]=d
        cols.append(col)
        dres3.append((val[3558]%p-b3)%p); dres29.append((val[29322]%p-b29)%p)
        val[w]=o
    forward(); ns['v']=val
    # build rows: for each currently-satisfied eq that some col touches, row = [col_j[i]], rhs=0
    touched=set()
    for col in cols: touched|=set(col)
    sat=[i for i in touched if baseroot[i]%p==0]  # keep these zero
    rows=[]; rhs=[]
    for i in sat:
        rows.append([cols[j].get(i,0) for j in range(len(W))]); rhs.append(0)
    rows.append(dres3); rhs.append((-b3)%p)
    rows.append(dres29); rhs.append((-b29)%p)
    x,rk=rref_solve(rows,rhs,len(W))
    if x is None:
        print(f"newton {newton}: INCONSISTENT (rank {rk}) - residues not reachable in witness null space"); break
    # apply mod-p step (balanced)
    for j,w in enumerate(W):
        d=x[j]%p
        if d> p//2: d-=p
        val[w]+=d
    forward(); ns['v']=val
    nb=sum(1 for i in range(len(lines)) if i not in CORE and eval(eqcode[i],ns)!=0)
    print(f"newton {newton}: applied; x3558%p={val[3558]%p==0}, x29322%p={val[29322]%p==0}, noncore_broken={nb}", flush=True)
# set quotients, final check
val[30317]=-(val[11150]//p); val[2936]=(537773*val[37758])//p
if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
forward(); ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); core={sorted(i for i in F if i in CORE)}")
print(f"noncore={sorted(i for i in F if i not in CORE)[:30]}")
print(f"S%p={val[35389]%p}, T%p={val[6671]%p}, L2/p%6672769={(val[25739]//p)%6672769}")
