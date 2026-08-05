#!/usr/bin/env python3
"""Test additivity of SPINE equation roots in the message bits, and build the full solve set."""
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
BITS=json.load(open('rootsys.json'))['BITS']
bits7=json.load(open('act7715.json'))['free7']
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
spine=[56,133,2071,2683,4386,7254,8073,11009,13660,15299,16622,17726,19066,19656,19712,20452,22093,25480,28090,28653,31061,31138,32894,34517,34892,35089,35299,38629]
scode=[compile(VAR.sub(r'v[\1]',inner_src(lines[i].rsplit('=',1)[0])),'<e>','eval') for i in spine]
ns={'__builtins__':{}}
def roots(setbits):
    val=val0[:]
    for b in setbits: val[b]=1
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    return [eval(c,ns) for c in scode]
r0=roots([])
print("spine roots at all-zero (which are nonzero):")
for k,i in enumerate(spine):
    if r0[k]!=0: print(f"  eq[{i}]: root=2^{r0[k].bit_length()-1 if r0[k]>0 else (abs(r0[k]).bit_length()-1)} ({'+' if r0[k]>0 else '-'})")
# additivity test on spine for several pairs
import itertools
print("\nspine additivity (violations per pair):")
viol_total=0
for Abit,Bbit in [(bits7[0],bits7[1]),(bits7[2],bits7[5]),(bits7[0],bits7[10]),(bits7[20],bits7[40])]:
    rA=roots([Abit]); rB=roots([Bbit]); rAB=roots([Abit,Bbit])
    v=[spine[k] for k in range(len(spine)) if rAB[k]!=rA[k]+rB[k]-r0[k]]
    print(f"  (x_{Abit},x_{Bbit}): {len(v)} violations {v[:6]}")
    viol_total+=len(v)
print(f"\nTOTAL spine additivity violations: {viol_total}")
