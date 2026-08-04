#!/usr/bin/env python3
"""Iterative integer solve for the feasible (1,1) config. Each round: F=failing eqs, accumulate
handles H (free inputs in F), Ffull=all eqs touching H, SNF-solve Ffull over H (roots->0), apply,
forward-eval. Repeat until 0 fail or no progress. Verify with checker."""
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
rootcode=[compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(L.rsplit('=',1)[0]))),'<e>','eval') for L in lines]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def snf_transform(Ain):
    m=len(Ain); n=len(Ain[0]); Awrk=[r[:] for r in Ain]
    U=[[1 if i==j else 0 for j in range(m)] for i in range(m)]
    V=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    def sr(i,j): Awrk[i],Awrk[j]=Awrk[j],Awrk[i]; U[i],U[j]=U[j],U[i]
    def sc(i,j):
        for r in Awrk: r[i],r[j]=r[j],r[i]
        for r in V: r[i],r[j]=r[j],r[i]
    def ar(i,j,f):
        for k in range(n): Awrk[i][k]+=f*Awrk[j][k]
        for k in range(m): U[i][k]+=f*U[j][k]
    def ac(i,j,f):
        for r in range(m): Awrk[r][i]+=f*Awrk[r][j]
        for r in range(n): V[r][i]+=f*V[r][j]
    for t in range(min(m,n)):
        while True:
            best=None;piv=None
            for i in range(t,m):
                for j in range(t,n):
                    if Awrk[i][j]!=0 and (best is None or abs(Awrk[i][j])<best): best=abs(Awrk[i][j]);piv=(i,j)
            if piv is None: return Awrk,U,V,t
            i0,j0=piv
            if i0!=t: sr(i0,t)
            if j0!=t: sc(j0,t)
            done=True
            for i in range(t+1,m):
                if Awrk[i][t]!=0:
                    ar(i,t,-(Awrk[i][t]//Awrk[t][t]))
                    if Awrk[i][t]!=0: done=False
            for j in range(t+1,n):
                if Awrk[t][j]!=0:
                    ac(j,t,-(Awrk[t][j]//Awrk[t][t]))
                    if Awrk[t][j]!=0: done=False
            if done:
                bad=False
                for i in range(t+1,m):
                    for j in range(t+1,n):
                        if Awrk[i][j]%Awrk[t][t]!=0: ar(t,i,1); bad=True; break
                    if bad: break
                if not bad: break
    return Awrk,U,V,min(m,n)
def solve_local(H, F):
    """Solve roots of F == 0 over handles H (integer). Returns True if applied a feasible step."""
    n=len(H); m=len(F)
    ns['v']=val
    base={i:eval(rootcode[i],ns) for i in F}
    Jac=[[0]*n for _ in range(m)]
    for j,h in enumerate(H):
        old=val[h]; val[h]=old+1; forward(); ns['v']=val
        for ri,i in enumerate(F): Jac[ri][j]=eval(rootcode[i],ns)-base[i]
        val[h]=old
    forward()
    D,U,V,rk=snf_transform([[Jac[i][j] for j in range(n)] for i in range(m)])
    bb=[-base[F[i]] for i in range(m)]
    c=[sum(U[i][k]*bb[k] for k in range(m)) for i in range(m)]
    yv=[0]*n
    for i in range(m):
        d=D[i][i] if i<n else 0
        if d==0:
            if c[i]!=0: return False,0
        else:
            if c[i]%d!=0: return False,0
            yv[i]=c[i]//d
    xint=[sum(V[i][j]*yv[j] for j in range(n)) for i in range(n)]
    for j,h in enumerate(H): val[h]+=xint[j]
    forward()
    return True, rk
forward(); ns['v']=val
H=set()
for it in range(10):
    ns['v']=val
    Ffail=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    print(f"iter {it}: {len(lines)-len(Ffail)}/{len(lines)} ({len(Ffail)} fail)", flush=True)
    if not Ffail:
        print("SOLVED!"); break
    # accumulate handles from failing eqs
    newH=set().union(*[eqvars[i]&hbase for i in Ffail])
    H |= newH
    Hlist=sorted(H)
    Ffull=[i for i in range(len(lines)) if eqvars[i] & H]
    print(f"  handles={len(Hlist)}, local system={len(Ffull)} eqs", flush=True)
    ok,rk=solve_local(Hlist, Ffull)
    if not ok:
        print(f"  infeasible at this scope (rank issue); stopping"); break
ns['v']=val
Ffail=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(Ffail)}/{len(lines)} ({len(Ffail)} fail): {Ffail[:20]}", flush=True)
if len(Ffail)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('iter_solved.json','w')); print("SAVED iter_solved.json")
