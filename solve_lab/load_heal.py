#!/usr/bin/env python3
"""From 39021 (core solved, 12 load-fails): set x_33462=CONST1, x_22152=CONST2, then Dixon-heal the
resulting breaks with handles OUTSIDE the two loads' forward cones (agent E: rank-9 consistent, 222
handles no feedback). Verify. If closes -> 39033 FULL SOLVE."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
C1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
C2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
s={int(k[2:]):v for k,v in json.load(open('best_agentA_39021.json')).items()}
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
val=[0]*NVARS
for k,v in s.items():
    if k<NVARS: val[k]=v
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp: ready[v]=True
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
# forward cones of x_33462, x_22152 (what they affect) -> exclude those vars' consumers
consumers=defaultdict(list)
for k,t in enumerate(order):
    for u in gates[definer[t]][2]: consumers[u].append(t)
def fwd_cone(root):
    seen=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        for t in consumers.get(x,()): st.append(t)
    return seen
loadcone=fwd_cone(33462)|fwd_cone(22152)|{33462,22152}
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
forward()
val[33462]=C1; val[22152]=C2
forward()
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"after set loads: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
def inv(a): return pow(a%p,p-2,p)
def gfp_piv(J):
    m=len(J); n=len(J[0]) if m else 0
    M=[[J[i][j]%p for j in range(n)] for i in range(m)]; rm=list(range(m)); pr=0; PR=[]; PC=[]
    for c in range(n):
        piv=None
        for i in range(pr,m):
            if M[i][c]%p: piv=i;break
        if piv is None: continue
        M[pr],M[piv]=M[piv],M[pr]; rm[pr],rm[piv]=rm[piv],rm[pr]
        iv=inv(M[pr][c]); M[pr]=[(x*iv)%p for x in M[pr]]
        for i in range(m):
            if i!=pr and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][k]-f*M[pr][k])%p for k in range(n)]
        PR.append(rm[pr]); PC.append(c); pr+=1
        if pr>=m: break
    return PR,PC
def matinv(M):
    r=len(M); Aug=[[M[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
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
def dixon(M,b,steps=16):
    r=len(M); Mi=matinv(M)
    if Mi is None: return None
    x=[0]*r; bb=b[:]; mod=1
    for _ in range(steps):
        bm=[bb[i]%p for i in range(r)]
        xi=[sum(Mi[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        nb=[]
        for i in range(r):
            ss=bb[i]-sum(M[i][k]*xi[k] for k in range(r))
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
Facc=set()
for rnd in range(15):
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if not F: print("SOLVED!"); break
    Facc|=set(F); Fl=sorted(Facc)
    H=sorted((set().union(*[eqvars[i]&freeinp for i in Fl])) - loadcone - {33462,22152})
    base=[eval(rc(i),ns) for i in Fl]
    Jac=[[0]*len(H) for _ in Fl]
    for j,h in enumerate(H):
        o=val[h]; val[h]=o+1; forward()
        for ri in range(len(Fl)): Jac[ri][j]=eval(rc(Fl[ri]),ns)-base[ri]
        val[h]=o
    forward()
    PR,PC=gfp_piv(Jac); r=len(PR)
    M=[[Jac[PR[i]][PC[j]] for j in range(r)] for i in range(r)]
    rhs=[-base[PR[i]] for i in range(r)]
    y=dixon(M,rhs)
    print(f"rnd {rnd}: {len(F)} now, {len(Fl)} accum, {len(H)} handles, rank {r}", flush=True)
    if y is None: print("  dixon fail"); break
    for j in range(r): val[H[PC[j]]]+=y[j]
    forward()
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {sorted(F)[:20]}")
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('FULL_SOLVED.json','w')); print("*** FULL SOLVED ***")
elif len(F)<12:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('load_heal_partial.json','w')); print(f"improved to {len(lines)-len(F)}")
