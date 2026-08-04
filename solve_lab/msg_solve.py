#!/usr/bin/env python3
"""Solve the message by linear algebra over ADDITIVE equation roots.
Additive rows = 16 large-root spine (|root|>2^50) + monitor rows (bit-affected). Non-additive
small-root spine (boolean control) excluded; handled by requiring x_7715=1 post-hoc.
System: sum_b m_b delta_{i,b} = -r0_i. Solve over Z (SNF), check boolean, forward-eval, verify."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
val0=[0]*NVARS; pinned=[False]*NVARS
for p in A:
    vs=atom_vars(p)
    if len(vs)==1:
        v=next(iter(vs)); c0=p.get((),0); c1=p.get((v,),0); c2=p.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val0[v]=(-c0)//c1; pinned[v]=True
rs=json.load(open('rootsys.json')); BITS=rs['BITS']; bidx={b:i for i,b in enumerate(BITS)}
gate_out=set(t for t,_,_ in gates)
for v in range(NVARS):
    if v not in gate_out: pinned[v]=True
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or pinned[v]: ready[v]=True
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
    return ast.unparse(node)
ns={'__builtins__':{}}
def roots_and_val(setbits, eqs):
    val=val0[:]
    for b in setbits: val[b]=1
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    codes=[compile(VAR.sub(r'v[\1]',inner_src(lines[i].rsplit('=',1)[0])),'<e>','eval') for i in eqs]
    return [eval(c,ns) for c in codes], val
spine_all=[56,133,2071,2683,4386,7254,8073,11009,13660,15299,16622,17726,19066,19656,19712,20452,22093,25480,28090,28653,31061,31138,32894,34517,34892,35089,35299,38629]
# baseline roots for spine
r0s,_=roots_and_val([], spine_all)
additive_spine=[spine_all[k] for k in range(len(spine_all)) if abs(r0s[k])>2**50]  # large-root
print(f"additive (large-root) spine: {len(additive_spine)} -> {additive_spine}", flush=True)
# delta for additive spine over all bits
codes_sp=[compile(VAR.sub(r'v[\1]',inner_src(lines[i].rsplit('=',1)[0])),'<e>','eval') for i in additive_spine]
def sp_roots(setbits):
    val=val0[:]
    for b in setbits: val[b]=1
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    return [eval(c,ns) for c in codes_sp]
r0=sp_roots([])
DELTA=[dict() for _ in additive_spine]
for i,b in enumerate(BITS):
    rb=sp_roots([b])
    for k in range(len(additive_spine)):
        d=rb[k]-r0[k]
        if d: DELTA[k][i]=d
    if i%60==0: print(f"  spine-delta {i}/256", flush=True)
# Build integer system rows: additive spine only (16 rows x 256). RHS = -r0
rows=[]; rhs=[]
for k in range(len(additive_spine)):
    row=[DELTA[k].get(j,0) for j in range(len(BITS))]
    rows.append(row); rhs.append(-r0[k])
# SNF solve
def snf_solve(Am,b):
    m=len(Am); nn=len(Am[0]) if Am else 0
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
            p=M[piv][piv]; done=True
            for i in range(m):
                if i!=piv and M[i][piv]!=0:
                    qq=M[i][piv]//p
                    for j in range(nn): M[i][j]-=qq*M[piv][j]
                    bb[i]-=qq*bb[piv]
                    if M[i][piv]!=0: done=False
            for j in range(nn):
                if j!=piv and M[piv][j]!=0:
                    qq=M[piv][j]//p
                    for i in range(m): M[i][j]-=qq*M[i][piv]
                    for i in range(nn): V[i][j]-=qq*V[i][piv]
                    if M[piv][j]!=0: done=False
            if done: break
    y=[0]*nn
    for i in range(min(m,nn)):
        d=M[i][i]
        if d==0: continue
        if bb[i]%d!=0: return None,('divfail',i,d,bb[i])
        y[i]=bb[i]//d
    for i in range(m):
        if all(z==0 for z in M[i]) and bb[i]!=0: return None,('inconsistent',i)
    return [sum(V[i][j]*y[j] for j in range(nn)) for i in range(nn)], None
sol,err=snf_solve(rows,rhs)
if sol is None:
    print(f"NO solution to additive-spine system: {err}", flush=True); sys.exit(0)
onbits=[BITS[j] for j in range(len(BITS)) if sol[j]!=0]
nonbool=[(BITS[j],sol[j]) for j in range(len(BITS)) if sol[j] not in (0,1)]
print(f"solution: {len(onbits)} bits set; non-boolean entries: {len(nonbool)}: {nonbool[:8]}", flush=True)
# forward-eval FULL with these bits (boolean rounding: set nonzero->1? first try exact)
setb=[BITS[j] for j in range(len(BITS)) if sol[j]==1]
allr, val = roots_and_val(setb, list(range(len(lines))))
fails=[i for i in range(len(lines)) if allr[i]!=0]
print(f"EXACT-bit (only sol==1) forward-eval: {len(lines)-len(fails)}/{len(lines)} fail={len(fails)}: {fails[:20]}", flush=True)
json.dump({'sol':[str(x) for x in sol],'onbits':onbits}, open('msg_sol.json','w'))
