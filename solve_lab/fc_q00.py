#!/usr/bin/env python3
"""Test quadrant (0,0): x_15298=0 trivializes the core. Run pure greedy forward-eval (no activators)
and forward_construct-style greedy for several override configs; report x_15298 and fail count."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates)
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
def run(override, label):
    val=[0]*NVARS; pinned=[False]*NVARS
    for pp in A:
        vs=atom_vars(pp)
        if len(vs)==1:
            v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
            if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
    freeinp=set(v for v in range(NVARS) if v not in gate_out)
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
    ns={'__builtins__':{}}; ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    cf=len([i for i in F if i in CORE])
    print(f"{label}: x_15298={val[15298]} x_7715={val[7715]} x_34554={val[34554]} | fails={len(F)} (core {cf}/20)")
    print(f"   noncore fails sample: {[i for i in F if i not in CORE][:12]}")
    return val,F
run({}, "no-override (natural, all free=0)")
run({30213:C2,22162:C1,24468:C1,18956:C2}, "constants-only (0,0)")
