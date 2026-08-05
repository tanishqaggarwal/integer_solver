#!/usr/bin/env python3
"""Robust iterative absorption. Activate + set forced constants + MUX data, then repair every
failing equation by adjusting a FREE INPUT that appears linearly in its root E (Gauss-Seidel
over Z). Free inputs feeding fewest equations preferred. Product gadgets zeroed by nulling a
free factor. Iterate to convergence; verify at equation level."""
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
for p in A:
    vs=atom_vars(p)
    if len(vs)==1:
        v=next(iter(vs)); c0=p.get((),0); c1=p.get((v,),0); c2=p.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
gate_out=set(t for t,_,_ in gates)
freeinp=set(v for v in range(NVARS) if v not in gate_out)
# pins: activation + forced constants + MUX data (deep free inputs)
override={ACT:1, 16742:C2, 12186:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
# topo (all free inputs are ready; targets computed)
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
rootnodes=[inner_src(L.rsplit('=',1)[0]) for L in lines]
rootcode=[compile(VAR.sub(r'v[\1]',ast.unparse(n)),'<e>','eval') for n in rootnodes]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
# free inputs per equation (candidate absorbers), and their global frequency
eq_free=[sorted((set(int(m) for m in VAR.findall(L)) & freeinp) - set(override)) for L in lines]
freq=defaultdict(int)
for s in eq_free:
    for v in s: freq[v]+=1
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def root(i):
    ns['v']=val; return eval(rootcode[i],ns)
# linear coefficient of free input v in equation i's root: finite diff
def coef(i,v):
    old=val[v]; base=root(i); val[v]=old+1; forward(); c=root(i)-base; val[v]=old; forward(); return c,base
forward()
for it in range(80):
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    if it%5==0 or len(F)<30: print(f"iter {it}: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
    if not F: print("SOLVED!"); break
    progress=False
    # sort free inputs by frequency (rare first) globally; repair each failing eq
    for i in F:
        base=root(i)
        if base==0: continue
        # try free inputs in this eq, rarest first, pick one with c|base
        chosen=None
        for v in sorted(eq_free[i], key=lambda z:freq[z]):
            c,_=coef(i,v)
            if c!=0 and base % c==0:
                chosen=(v,base//c); break
        if chosen:
            v,delta=chosen; val[v]-=delta; forward(); progress=True
    if not progress:
        print("  no linear absorber found for remaining fails; break"); break
ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:30]}", flush=True)
if len(F)<26:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('absorb_solved.json','w'))
    print("saved absorb_solved.json")
