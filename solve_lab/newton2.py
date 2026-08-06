#!/usr/bin/env python3
"""Simultaneous absorber solve WITH forced constants pinned. Pin activation+forced-constants+
MUX data; solve remaining free-input handles (excluding bits) over equation ROOTS via SNF."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
ACT=int(sys.argv[1]) if len(sys.argv)>1 else 22106
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
override={ACT:1, 16742:C2, 12186:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
handles_all=sorted(freeinp - boolbits - set(override)); hset=set(handles_all)
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
rootcode=[compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(L.rsplit('=',1)[0]))),'<e>','eval') for L in lines]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def snf_solve(Am,b):
    m=len(Am); nn=len(Am[0]) if Am else 0
    if m==0: return [0]*nn,None
    M=[r[:] for r in Am]; bb=b[:]
    V=[[1 if i==j else 0 for j in range(nn)] for i in range(nn)]
    for piv in range(min(m,nn)):
        while True:
            best=None;found=None
            for i in range(piv,m):
                for j in range(piv,nn):
                    if M[i][j]!=0 and (best is None or abs(M[i][j])<best): best=abs(M[i][j]);found=(i,j)
            if found is None: break
            i0,j0=found
            if i0!=piv: M[piv],M[i0]=M[i0],M[piv]; bb[piv],bb[i0]=bb[i0],bb[piv]
            if j0!=piv:
                for r in M: r[piv],r[j0]=r[j0],r[piv]
                for r in V: r[piv],r[j0]=r[j0],r[piv]
            ppv=M[piv][piv]; done=True
            for i in range(m):
                if i!=piv and M[i][piv]!=0:
                    qq=M[i][piv]//ppv
                    for j in range(nn): M[i][j]-=qq*M[piv][j]
                    bb[i]-=qq*bb[piv]
                    if M[i][piv]!=0: done=False
            for j in range(nn):
                if j!=piv and M[piv][j]!=0:
                    qq=M[piv][j]//ppv
                    for i in range(m): M[i][j]-=qq*M[i][piv]
                    for i in range(nn): V[i][j]-=qq*V[i][piv]
                    if M[piv][j]!=0: done=False
            if done: break
    y=[0]*nn
    for i in range(min(m,nn)):
        d=M[i][i]
        if d==0: continue
        if bb[i]%d!=0: return None,('divfail',i,d)
        y[i]=bb[i]//d
    for i in range(m):
        if all(z==0 for z in M[i]) and bb[i]!=0: return None,('inconsistent',i)
    return [sum(V[i][j]*y[j] for j in range(nn)) for i in range(nn)],None
forward()
for it in range(15):
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
    if not F: print("SOLVED!"); break
    H=sorted(set().union(*[eqvars[i]&hset for i in F]))
    base=[eval(rootcode[i],ns) for i in F]
    J=[[0]*len(H) for _ in F]
    for j,h in enumerate(H):
        old=val[h]; val[h]=old+1; forward(); ns['v']=val
        for ri,i in enumerate(F): J[ri][j]=eval(rootcode[i],ns)-base[ri]
        val[h]=old
    forward()
    d,err=snf_solve(J,[-x for x in base])
    if d is None:
        print(f"  linear solve failed: {err} at eq[{F[err[1]]}] (nhandles={len(H)})"); break
    for j,h in enumerate(H): val[h]+=d[j]
    forward()
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:30]}")
if len(F)<26:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('newton2_solved.json','w')); print("SAVED")
