#!/usr/bin/env python3
"""Solve the 256-bit message by LINEAR ALGEBRA over equation roots.
Each equation's root E (strip square/const) is linear in the bits (additivity verified).
Want E_i(m)=0 for all monitor eqs. Build delta_{i,b}=E_i(e_b)-E_i(0); solve
sum_b m_b delta_{i,b} = -E_i(0) over Z, check boolean, forward-eval, verify."""
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
bits7=json.load(open('act7715.json'))['free7']
# x_34554 cone bits
gate_defs=defaultdict(list)
for t,rhs,vids in gates: gate_defs[t].append(rhs)
summ=json.load(open('atoms/summary.json')); inputs=set(summ['inputs'])
def cone(v, seen):
    if v in seen: return
    seen.add(v)
    for rhs in gate_defs.get(v,[]):
        for m in re.findall(r'x_(\d+)',rhs): cone(int(m),seen)
sys.setrecursionlimit(1000000)
c34=set(); cone(34554,c34)
bits34=[v for v in c34 if v in inputs or v not in gate_defs]
BITS=sorted(set(bits7)|set(bits34))
print(f"total message bits: {len(BITS)} (x_7715 cone {len(bits7)}, x_34554 cone {len(bits34)})", flush=True)
gate_out=set(t for t,_,_ in gates)
for v in range(NVARS):
    if v not in gate_out: pinned[v]=True
# topo
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
# root (inner E) for each equation
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
bitset=set(BITS)
mon=[i for i,L in enumerate(lines) if set(int(m) for m in VAR.findall(L)) & bitset]
rootcode=[compile(VAR.sub(r'v[\1]',inner_src(lines[i].rsplit('=',1)[0])),'<e>','eval') for i in mon]
print(f"monitor eqs: {len(mon)}", flush=True)
ns={'__builtins__':{}}
def roots(setbits):
    val=val0[:]
    for b in setbits: val[b]=1
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    return [eval(c,ns) for c in rootcode], val
r0,_=roots([])
# spine = monitor eqs with nonzero root at all-zero
spine=[k for k in range(len(mon)) if r0[k]!=0]
print(f"spine (nonzero root at all-zero): {len(spine)} -> {[mon[k] for k in spine]}", flush=True)
# delta matrix (only rows where some bit matters); compute per-bit
delta=[dict() for _ in range(len(mon))]  # delta[k][bidx]
for bidx,b in enumerate(BITS):
    rb,_=roots([b])
    for k in range(len(mon)):
        d=rb[k]-r0[k]
        if d!=0: delta[k][bidx]=d
    if bidx%50==0: print(f"  signatures {bidx}/{len(BITS)}", flush=True)
# rows that any bit affects OR that are spine
rows=[k for k in range(len(mon)) if delta[k] or r0[k]!=0]
print(f"active rows (bit-affected or spine): {len(rows)}", flush=True)
json.dump({'BITS':BITS,'mon':mon,'rows':rows,'r0':[str(r0[k]) for k in range(len(mon))],
           'delta':[{str(j):str(v) for j,v in delta[k].items()} for k in rows]}, open('rootsys.json','w'))
print("saved rootsys.json", flush=True)
