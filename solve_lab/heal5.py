#!/usr/bin/env python3
"""From best + minimal residue fix (S=T=0 mod p) + quotient handles: Dixon-solve the broken NON-CORE
wiring equations using handles disjoint from the residue-critical set (x_24908 cone, controls,
quotients). Keeps S,T=0 mod p while re-healing wiring."""
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
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
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
def freecone(root):
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    return leaves
def inv_modp(a): return pow(a%p,p-2,p)
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
        iv=inv_modp(M[pr][c]); M[pr]=[(x*iv)%p for x in M[pr]]
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
        Aug[c],Aug[piv]=Aug[piv],Aug[c]; iv=inv_modp(Aug[c][c]); Aug[c]=[(x*iv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p!=0:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(M,b,steps=8):
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
# ripple cost: #equations a handle appears in outside broken+core (computed lazily on full eq set)
byh=defaultdict(set)
for _i,_vs in enumerate(eqvars):
    for _v in _vs&freeinp: byh[_v].add(_i)
CORE_=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward()
val[14853]=val[14853]-(val[29322]%p); forward()
val[16742]=val[16742]+(val[3558]%p); forward()
val[30317]=-(val[11150]//p); val[2936]=(537773*val[37758])//p; val[5146]=val[25739]//(6672769*p); forward()
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
crit={14853,12186,16742,30317,2936,5146}|freecone(24908)
ns['v']=val
Fall=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
Fnc=[i for i in Fall if i not in CORE]
print(f"start: {len(lines)-len(Fall)}/{len(lines)}; non-core broken={len(Fnc)}: {Fnc}")
for rnd in range(15):
    val[30317]=-(val[11150]//p); val[2936]=(537773*val[37758])//p; val[5146]=val[25739]//(6672769*p); forward()
    ns['v']=val
    Fnc=[i for i in range(len(lines)) if i not in CORE and eval(eqcode[i],ns)!=0]
    if not Fnc: print(f"round {rnd}: all non-core healed"); break
    Fncset=set(Fnc)
    Hset=set().union(*[eqvars[i]&freeinp for i in Fnc]) - crit
    def ripple(h): return len(byh[h]-Fncset-CORE_)
    H=sorted(Hset, key=ripple)
    print(f"round {rnd}: {len(Fnc)} broken, {len(H)} handles (min ripple {ripple(H[0]) if H else 0})", flush=True)
    base=[eval(rootcode_of(i),ns) for i in Fnc]
    Jac=[[0]*len(H) for _ in range(len(Fnc))]
    for j,h in enumerate(H):
        o=val[h]; val[h]=o+1; forward(); ns['v']=val
        for ri in range(len(Fnc)): Jac[ri][j]=eval(rootcode_of(Fnc[ri]),ns)-base[ri]
        val[h]=o
    forward(); ns['v']=val
    pr,pc=gfp_pivots(Jac); r=len(pr)
    M=[[Jac[pr[i]][pc[j]] for j in range(r)] for i in range(r)]
    rhs=[-base[pr[i]] for i in range(r)]
    y=dixon(M,rhs)
    if y is None: print(f"  dixon failed (rank {r}/{len(Fnc)})"); break
    for j in range(r): val[H[pc[j]]]+=y[j]
    forward(); ns['v']=val
    print(f"  applied rank {r}; S%p={val[35389]%p==0}, T%p={val[6671]%p==0}")
ns['v']=val
Fall=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(Fall)}/{len(lines)} ({len(Fall)} fail)")
print(f"core fail: {sorted(i for i in Fall if i in CORE)}")
print(f"non-core fail: {sorted(i for i in Fall if i not in CORE)}")
print(f"S%p={val[35389]%p}, T%p={val[6671]%p}, L2/p%6672769={(val[25739]//p)%6672769}")
if all(i in CORE for i in Fall):
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('heal_state.json','w')); print("SAVED heal_state.json (all non-core satisfied)")
