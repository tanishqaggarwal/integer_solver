#!/usr/bin/env python3
"""Solve the message via INHOMOGENEOUS GF(p) conditions. From build_twist (x_22106 active +
forced constants + MUX data), the 27 failing roots must become 0 (mod p first, then Z). Extra
x_7715-cone bits (keep x_7715=1) contribute additively. Solve sum_b m_b*delta_b == -root0 (mod p)
over GF(p), search for a boolean solution, then let data absorb the Z-quotient. Verify."""
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
val0=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val0[v]=(-c0)//c1; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
bits7=json.load(open('act7715.json'))['free7']
override={22106:1, 16742:C2, 12186:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val0[v]=x
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
ns={'__builtins__':{}}
def run(extra):
    val=val0[:]
    for b in extra: val[b]=1
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    return val
val=run([]); ns['v']=val
FAIL=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"build_twist failing: {len(FAIL)}", flush=True)
rootcode={i:compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval') for i in FAIL}
def roots(val):
    ns['v']=val
    return {i:eval(rootcode[i],ns) for i in FAIL}
r0=roots(val)
# candidate extra bits: x_7715 cone bits not already set
extras=[b for b in bits7 if b!=22106]
# mod-p signature of each extra bit on each failing root
print(f"computing mod-p signatures for {len(extras)} extra bits...", flush=True)
SIG={}
for bi,b in enumerate(extras):
    v=run([b]); rb=roots(v)
    SIG[b]={i:(rb[i]-r0[i])%p for i in FAIL if (rb[i]-r0[i])%p!=0}
    if bi%40==0: print(f"  {bi}/{len(extras)}", flush=True)
# Build GF(p) system: sum_b m_b*SIG[b][i] == (-r0[i]) mod p, for each i in FAIL
rowsFAIL=[i for i in FAIL]
target={i:(-r0[i])%p for i in FAIL}
# matrix cols=extras
def inv(a): return pow(a%p,p-2,p)
# Gaussian elim over GF(p): rows=FAIL eqs, cols=extras; augmented with target
M=[]
for i in rowsFAIL:
    row=[SIG[b].get(i,0) for b in extras]+[target[i]]
    M.append(row)
ne=len(extras)
# eliminate
prow=0; pivots=[]
for c in range(ne):
    piv=None
    for r in range(prow,len(M)):
        if M[r][c]%p!=0: piv=r;break
    if piv is None: continue
    M[prow],M[piv]=M[piv],M[prow]
    iv=inv(M[prow][c]); M[prow]=[(x*iv)%p for x in M[prow]]
    for r in range(len(M)):
        if r!=prow and M[r][c]%p!=0:
            f=M[r][c]; M[r]=[(M[r][j]-f*M[prow][j])%p for j in range(ne+1)]
    pivots.append(c); prow+=1
    if prow>=len(M): break
# consistency: any row 0..0 | nonzero?
incons=any(all(M[r][j]%p==0 for j in range(ne)) and M[r][ne]%p!=0 for r in range(len(M)))
print(f"GF(p) system: {len(rowsFAIL)} conditions, {ne} bit-vars, rank={len(pivots)}, consistent={not incons}", flush=True)
if incons:
    print("INCONSISTENT mod p -> no message (even non-boolean) fixes it with this activation.")
    sys.exit(0)
# particular solution: free vars=0
sol={c:0 for c in range(ne)}
pivset=set(pivots)
for r,c in enumerate(pivots):
    sol[c]=M[r][ne]%p
# check boolean-ness of particular solution
nonbool=[extras[c] for c in range(ne) if sol[c] not in (0,1)]
onbits=[extras[c] for c in range(ne) if sol[c]==1]
print(f"particular GF(p) solution: {len(onbits)} bits=1, {len(nonbool)} non-boolean", flush=True)
json.dump({'pivots':pivots,'sol':{str(extras[c]):str(sol[c]) for c in range(ne)},
           'kerdim':ne-len(pivots),'onbits':onbits,'nonbool_count':len(nonbool)}, open('p_msg.json','w'))
# if particular solution is boolean, test it
if not nonbool:
    v=run(onbits); ns['v']=v
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    print(f"BOOLEAN message found! forward-eval: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)", flush=True)
    if len(F)<26: json.dump({f"x_{i}":v[i] for i in range(NVARS)}, open('pmsg_solved.json','w'))
