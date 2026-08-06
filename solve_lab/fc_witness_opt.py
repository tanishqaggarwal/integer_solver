#!/usr/bin/env python3
"""Optimized witness null-space solve targeting S=x_35389, T=x_6671 residues mod p.
Partial-forward (only downstream gates) for speed. Expanded 2-hop witness. Reports consistency."""
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
# ancestors (free-input) and gate topo-position index
posof={t:k for k,t in enumerate(order)}
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
ns={'__builtins__':{}}; ns['v']=val
def forward():
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
# gate consumers for downstream cone
consumers=defaultdict(list)
for k,t in enumerate(order):
    for u in gates[definer[t]][2]:
        consumers[u].append(k)
def downstream_ks(w):
    # topo positions of gates affected by free input w
    seen=set(); st=[w]; res=[]
    # BFS over consumers, collect gate order positions
    affected=set()
    stack=[w]
    seenv=set()
    while stack:
        x=stack.pop()
        if x in seenv: continue
        seenv.add(x)
        for k in consumers.get(x,()):
            t=order[k]
            if k not in affected:
                affected.add(k); stack.append(t)
    return sorted(affected)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward()
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
# 2-hop witness expansion
W0=set(v for v in freeinp if best.get(v,0)!=0 and v not in override)
byh=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]&freeinp: byh[v].add(i)
extra=set()
for i in range(len(lines)):
    if eqvars[i]&{14853,12186,16742}: extra|=(eqvars[i]&freeinp)
def fcone(root):
    seen=set(); lv=set(); st=[root]
    gd={t:vids for t,rhs,vids in gates}
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gd:
            for u in gd[x]: st.append(u)
        elif x in freeinp: lv.add(x)
    return lv
extra|=fcone(24908)
hop1=set()
for h in (extra|{14853,12186,16742}): hop1|=byh[h]
for i in hop1: extra|=(eqvars[i]&freeinp)
W=sorted((W0|extra)-set(override))
print(f"witness inputs: {len(W)}", flush=True)
# precompute downstream ks and affected equations per W
dks={}; deq={}
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
for w in W:
    ks=downstream_ks(w)
    dks[w]=ks
    aff=set(eqbyvar.get(w,set()))
    for k in ks: aff|=eqbyvar.get(order[k],set())
    deq[w]=sorted(aff)
posS=posof[35389]; posT=posof[6671]
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
for newton in range(10):
    forward()
    bS=val[35389]%p; bT=val[6671]%p
    if bS==0 and bT==0: print(f"newton {newton}: S,T ZERO mod p"); break
    baseroot={i:eval(rootcode_of(i),ns)%p for i in range(len(lines))}
    cols=[]; dS=[]; dT=[]
    for wi,w in enumerate(W):
        o=val[w]; val[w]=o+1
        for k in dks[w]: val[order[k]]=eval(gcode[k],ns)
        col={}
        for i in deq[w]:
            d=(eval(rootcode_of(i),ns)-baseroot[i])%p
            if d: col[i]=d
        cols.append(col); dS.append((val[35389]%p-bS)%p); dT.append((val[6671]%p-bT)%p)
        val[w]=o
        for k in dks[w]: val[order[k]]=eval(gcode[k],ns)
    touched=set()
    for col in cols: touched|=set(col)
    sat=[i for i in touched if baseroot[i]==0]
    rows=[]; rhs=[]
    for i in sat:
        rows.append([cols[j].get(i,0) for j in range(len(W))]); rhs.append(0)
    rows.append(dS); rhs.append((-bS)%p)
    rows.append(dT); rhs.append((-bT)%p)
    x,rk=rref_solve(rows,rhs,len(W))
    print(f"newton {newton}: |sat|={len(sat)}, rank={rk}, {'CONSISTENT' if x is not None else 'INCONSISTENT'}", flush=True)
    if x is None: break
    for j,w in enumerate(W):
        d=x[j]%p
        if d>p//2: d-=p
        val[w]+=d
    forward()
    nb=sum(1 for i in range(len(lines)) if i not in CORE and eval(eqcode[i],ns)!=0)
    print(f"  applied; S%p={val[35389]%p==0}, T%p={val[6671]%p==0}, noncore_broken={nb}", flush=True)
val[30317]=-(val[11150]//p); val[2936]=(537773*val[37758])//p
if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
forward()
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); core={sorted(i for i in F if i in CORE)}")
print(f"S%p={val[35389]%p}, T%p={val[6671]%p}")
if len([i for i in F if i not in CORE])==0 and len(F)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('witness_solved.json','w')); print("SAVED witness_solved.json")
