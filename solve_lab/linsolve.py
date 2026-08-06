#!/usr/bin/env python3
"""Fix control bits (activator=1 => x_7715=1,x_34554=0). Then the failing system should be
LINEAR in the free data handles. Build Jacobian by finite differences over forward-eval,
solve J*h = -b over Z (with the huge pins as unknowns too), verify."""
import json, re, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
ACT=47
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
val[ACT]=1; pinned[ACT]=True
# handles = data free inputs to solve for
handles=json.load(open('data_handles.json'))
for h in [24468,18956]:
    if h not in handles: handles.append(h)
hidx={h:i for i,h in enumerate(handles)}
for h in handles: pinned[h]=True   # treat as pinned (we set them)
# topo orientation
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
gcode={t:compile(VAR.sub(r'v[\1]',gates[definer[t]][1]),'<r>','eval') for t in order}
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
# only the failing eqs matter, but re-derive fails at current handle=0
ecode={i:compile(VAR.sub(r'v[\1]',lines[i].rsplit('=',1)[0]),'<e>','eval') for i in range(len(lines))}
ns={'__builtins__':{}}
def forward(hv):
    for h,i in hidx.items(): val[h]=hv[i]
    ns['v']=val
    for t in order: val[t]=eval(gcode[t],ns)
def all_resid():
    ns['v']=val
    return [eval(ecode[i],ns) for i in range(len(lines))]
# baseline with all handles=0
h0=[0]*len(handles)
forward(h0); base=all_resid()
FAIL=[i for i,r in enumerate(base) if r!=0]
print(f"failing at handles=0: {len(FAIL)}", flush=True)
b=[base[i] for i in FAIL]
# Jacobian: perturb each handle by 1
J=[[0]*len(handles) for _ in FAIL]
for j,h in enumerate(handles):
    hv=h0[:]; hv[j]=1; forward(hv); r=all_resid()
    for r_i,i in enumerate(FAIL): J[r_i][j]=r[i]-base[i]
    # linearity check on first handle
forward(h0)
# check linearity: handle j=0 at value 2
hv=h0[:]; hv[0]=2; forward(hv); r2=all_resid()
lin_ok=all((r2[i]-base[i])==2*J[k][0] for k,i in enumerate(FAIL))
print(f"linearity check (handle0 @2 == 2*Jacobian): {lin_ok}", flush=True)
json.dump({'FAIL':FAIL,'b':[str(x) for x in b],'J':[[str(x) for x in row] for row in J],'handles':handles}, open('linsys.json','w'))
print(f"built linear system: {len(FAIL)} eqs x {len(handles)} unknowns", flush=True)

# --- per-handle linearity test ---
print("\n=== per-handle linearity (which handles are nonlinear) ===", flush=True)
forward(h0); base=all_resid()
nonlin=[]; lin=[]
for j,h in enumerate(handles):
    hv=h0[:]; hv[j]=1; forward(hv); r1=all_resid()
    hv[j]=2; forward(hv); r2=all_resid()
    d1=[r1[i]-base[i] for i in FAIL]; d2=[r2[i]-base[i] for i in FAIL]
    islin=all(d2[k]==2*d1[k] for k in range(len(FAIL)))
    (lin if islin else nonlin).append(h)
forward(h0)
print(f"linear handles: {len(lin)}; nonlinear handles: {len(nonlin)}")
print(f"nonlinear: {['x_%d'%h for h in nonlin]}")
json.dump({'lin':lin,'nonlin':nonlin}, open('handle_class.json','w'))
