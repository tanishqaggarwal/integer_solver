#!/usr/bin/env python3
"""Agent A's proposed fix: from 39013, zero x_3558 by moving x_24908 (via an additive free-input knob
in its cone) to equal x_16742, and zero x_29322 by moving x_14853 = x_12186 — keeping x_16742 and
x_12186 (the load-partners of x_33462/x_22152) at baseline so the 12 loads stay intact. Then set
quotient handles + verify."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={t:vids for t,rhs,vids in gates}
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
val=[0]*NVARS
for k,v in best.items():
    if k<NVARS: val[k]=v
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); pinned=[False]*NVARS
for v in freeinp: pinned[v]=True
ready=[False]*NVARS
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
def freecone(root):
    seen=set(); lv=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x]: st.append(u)
        elif x in freeinp: lv.add(x)
    return lv
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{},'v':val}
def forward():
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
forward()
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"39013 baseline: {len(F0)} fail; x_3558%p!=0:{val[3558]%p!=0}, x_29322%p!=0:{val[29322]%p!=0}")
# find an additive knob for x_24908: a free input in its cone with d(x_24908)=1
xc=freecone(24908)
knob=None
b24908=val[24908]
for v in sorted(xc):
    o=val[v]; val[v]=o+1; forward()
    d=(val[24908]-b24908)
    val[v]=o
    if d==1: knob=v; break
    forward()
forward()
print(f"x_24908 additive-1 knob: x_{knob}")
if knob is None:
    # fallback: any knob with small coeff, or use x_16742 side
    print("no coeff-1 knob; aborting")
    sys.exit()
# zero x_3558 = x_24908 - x_16742 -> set x_24908 = x_16742 via knob (mod p, keep x_16742 baseline)
target=val[16742]
# need x_24908 -> target. current x_24908=b24908. delta = (target - b24908). knob has coeff1 so knob += delta%p for mod-p zero
val[knob]= val[knob] + (target - val[24908])%p
forward()
print(f"after x_24908 knob: x_3558%p={val[3558]%p}, x_24908%p match x_16742%p: {val[24908]%p==val[16742]%p}")
# zero x_29322 = x_14853 - x_12186 -> set x_14853 = x_12186 (mod p), keep x_12186 baseline
val[14853]= val[14853] - (val[29322]%p)
forward()
print(f"after x_14853: x_29322%p={val[29322]%p}, S%p={val[35389]%p}, T%p={val[6671]%p}")
# set quotient handles
if val[11150]%p==0: val[30317]=-(val[11150]//p)
if val[37758]%p==0: val[2936]=(537773*val[37758])//p
if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
forward()
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); core={len([i for i in F if i in CORE])}, noncore={sorted(i for i in F if i not in CORE)[:15]}")
if len(F)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('FULL_SOLVED.json','w')); print("*** FULL SOLVED ***")
elif len(F)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('agentAfix_partial.json','w')); print(f"saved partial {len(lines)-len(F)}")
