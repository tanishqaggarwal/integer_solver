#!/usr/bin/env python3
"""Two-phase: from wire1_m2fixed (core solved via wire=1, 27 noncore fail), Newton over
(wire-as-ONE-variable V + core quotients x_30317,x_5146,x_2936 + activated handles) to zero the 20
core AND 27 noncore jointly. Wire moves to heal the 13 unpackings; quotients co-adjust so core stays
0 (bilinear handled by including both in the linear system). Iterate. Verify."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977; MOD=6672769
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
r0,_=find2(26064); wsign={v:find2(v)[1] for v in range(NVARS) if find2(v)[0]==r0}
wire=set(wsign)
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={t:vids for t,rhs,vids in gates}
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
s0={int(k[2:]):v for k,v in json.load(open('wire1_m2fixed.json')).items()}
val=[0]*NVARS
for k,v in s0.items():
    if k<NVARS: val[k]=v
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp or v in wire: ready[v]=True
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
def rc(i):
    if i not in _rc: _rc[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rc[i]
ns={'__builtins__':{},'v':val}
def forward():
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def freecone(root):
    seen=set(); lv=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x]: st.append(u)
        elif x in freeinp: lv.add(x)
    return lv
loadcone=freecone(33462)|freecone(22152)|{33462,22152,16742,12186,19083}
forward()
def inv(a): return pow(a%p,p-2,p)
def gfp(J):
    m=len(J); n=len(J[0]) if m else 0; Mx=[[J[i][j]%p for j in range(n)] for i in range(m)]; rm=list(range(m)); pr=0; PR=[];PC=[]
    for c in range(n):
        piv=None
        for i in range(pr,m):
            if Mx[i][c]%p: piv=i;break
        if piv is None: continue
        Mx[pr],Mx[piv]=Mx[piv],Mx[pr]; rm[pr],rm[piv]=rm[piv],rm[pr]; iv=inv(Mx[pr][c]); Mx[pr]=[(x*iv)%p for x in Mx[pr]]
        for i in range(m):
            if i!=pr and Mx[i][c]%p:
                f=Mx[i][c]; Mx[i]=[(Mx[i][k]-f*Mx[pr][k])%p for k in range(n)]
        PR.append(rm[pr]); PC.append(c); pr+=1
        if pr>=m: break
    return PR,PC
def matinv(Mx):
    r=len(Mx); Aug=[[Mx[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv=None
        for i in range(c,r):
            if Aug[i][c]%p: piv=i;break
        if piv is None: return None
        Aug[c],Aug[piv]=Aug[piv],Aug[c]; iv=inv(Aug[c][c]); Aug[c]=[(x*iv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(Mx,b,steps=14):
    r=len(Mx); Mi=matinv(Mx)
    if Mi is None: return None
    x=[0]*r; bb=b[:]; mod=1
    for _ in range(steps):
        bm=[bb[i]%p for i in range(r)]; xi=[sum(Mi[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        nb=[]
        for i in range(r):
            ss=bb[i]-sum(Mx[i][k]*xi[k] for k in range(r))
            if ss%p: return None
            nb.append(ss//p)
        bb=nb; mod*=p
        if all(z==0 for z in bb): break
    half=mod//2; y=[]
    for xi in x:
        xi%=mod
        if xi>half: xi-=mod
        y.append(xi)
    return y
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
def perturb_wire(delta):
    for v,s in wsign.items(): val[v]+=s*delta
for it in range(6):
    forward()
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if not F: print("*** SOLVED ***"); break
    # equations to zero: all currently-failing + all core (keep solved)
    E=sorted(set(F)|CORE)
    # handles: WIRE-as-one (special index -1) + core quotients + activated free inputs in E cones minus loadcone
    H=[-1, 30317, 5146, 2936]
    extra=sorted((set().union(*[eqvars[i]&freeinp for i in E])) - loadcone - wire - set(H))
    H=H+extra
    base=[eval(rc(i),ns) for i in E]
    Jac=[[0]*len(H) for _ in E]
    for j,h in enumerate(H):
        if h==-1:
            perturb_wire(1); forward()
            for ri in range(len(E)): Jac[ri][j]=eval(rc(E[ri]),ns)-base[ri]
            perturb_wire(-1); forward()
        else:
            o=val[h]; val[h]=o+1; forward()
            for ri in range(len(E)): Jac[ri][j]=eval(rc(E[ri]),ns)-base[ri]
            val[h]=o
    forward()
    PR,PC=gfp(Jac); r=len(PR)
    Mx=[[Jac[PR[i]][PC[j]] for j in range(r)] for i in range(r)]; rhs=[-base[PR[i]] for i in range(r)]
    y=dixon(Mx,rhs)
    print(f"it {it}: {len(F)} noncore-fail, |E|={len(E)}, {len(H)} handles, rank {r}", flush=True)
    if y is None: print("  dixon fail"); break
    for j in range(r):
        h=H[PC[j]]
        if h==-1: perturb_wire(y[j])
        else: val[h]+=y[j]
    forward()
forward()
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {sorted(F)[:20]}")
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('FULL_SOLVED.json','w')); print("*** FULL SOLVED ***")
elif len(F)<27:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('twophase_partial.json','w')); print(f"improved to {len(lines)-len(F)}")
