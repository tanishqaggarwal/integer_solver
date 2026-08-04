#!/usr/bin/env python3
"""Test SNF integer-feasibility of the failing system across quadrant routings.
Quadrants: (1,0)->x_34606 routes x_16742,x_12186; (0,1)->x_5647 routes x_24908,x_14853;
(1,1)->x_15298 routes x_30213,x_22162. For each, activate the needed control cone(s), pin
forced constants + routed data, build 27x(handles) integer Jacobian, run SNF feasibility."""
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
val0=[0]*NVARS; pin0=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pin0[v]: val0[v]=(-c0)//c1; pin0[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
boolbits=set(json.load(open('boolbits.json'))['boolvars'])
act7=json.load(open('act7715.json'))['free7'][0]  # any x_7715 activator... use a good one
act7=22106
act34=json.load(open('act34554.json'))[0]
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
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
def snf_transform(Ain):
    m=len(Ain); n=len(Ain[0]); Awrk=[r[:] for r in Ain]
    U=[[1 if i==j else 0 for j in range(m)] for i in range(m)]
    def sr(i,j): Awrk[i],Awrk[j]=Awrk[j],Awrk[i]; U[i],U[j]=U[j],U[i]
    def sc(i,j):
        for r in Awrk: r[i],r[j]=r[j],r[i]
    def ar(i,j,f):
        for k in range(n): Awrk[i][k]+=f*Awrk[j][k]
        for k in range(m): U[i][k]+=f*U[j][k]
    def ac(i,j,f):
        for r in range(m): Awrk[r][i]+=f*Awrk[r][j]
    for t in range(min(m,n)):
        while True:
            best=None;piv=None
            for i in range(t,m):
                for j in range(t,n):
                    if Awrk[i][j]!=0 and (best is None or abs(Awrk[i][j])<best): best=abs(Awrk[i][j]);piv=(i,j)
            if piv is None: return Awrk,U,t
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
    return Awrk,U,min(m,n)
def test(name, override, acts):
    val=val0[:]; pinned=pin0[:]
    for a in acts: val[a]=1; pinned[a]=True
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
    gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
    ns={'__builtins__':{}}
    def fwd():
        ns['v']=val
        for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    fwd(); ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    hset=freeinp - boolbits - set(override) - set(acts)
    if not F:
        print(f"  {name}: 0 failing (already solved?!)"); return
    rc={i:compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval') for i in F}
    H=sorted(set().union(*[eqvars[i]&hset for i in F]))
    base={i:eval(rc[i],ns) for i in F}
    Jac=[[0]*len(H) for _ in F]
    for j,h in enumerate(H):
        old=val[h]; val[h]=old+1; fwd(); ns['v']=val
        for ri,i in enumerate(F): Jac[ri][j]=eval(rc[i],ns)-base[i]
        val[h]=old
    fwd()
    m=len(F); n=len(H)
    D,U,rk=snf_transform([[Jac[i][j] for j in range(n)] for i in range(m)])
    bb=[-base[F[i]] for i in range(m)]
    c=[sum(U[i][k]*bb[k] for k in range(m)) for i in range(m)]
    feas=True; obstr=0
    for i in range(m):
        d=D[i][i] if i<n else 0
        if d==0:
            if c[i]!=0: feas=False; obstr+=1
        elif c[i]%d!=0: feas=False; obstr+=1
    return feas,obstr,len(F)

act7s=json.load(open('act7715.json'))['free7'][:8]
act34s=json.load(open('act34554.json'))[:8]
best=(99,None)
print("scan (1,1):", flush=True)
for a7 in act7s:
    for a34 in act34s:
        feas,obstr,nf=test("s",{30213:C2,22162:C1,24468:C1,18956:C2},[a7,a34])
        if obstr<best[0]: best=(obstr,(a7,a34)); print(f"  a7={a7} a34={a34}: obstr={obstr} feasible={feas}", flush=True)
        if feas:
            print(f"  FEASIBLE at a7={a7} a34={a34}!", flush=True)
print("best:", best, flush=True)
