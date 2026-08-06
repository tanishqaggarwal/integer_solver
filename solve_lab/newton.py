#!/usr/bin/env python3
"""Newton-over-integers: fix activation bits, solve data/absorber free-inputs simultaneously.
From the activated state, iterate: forward-eval -> failing eqs -> Jacobian wrt free-input
handles -> SNF integer solve J*d=-r -> apply -> repeat."""
import json, re, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
ACT=int(sys.argv[1]) if len(sys.argv)>1 else 22106
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
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
boolbits=set(json.load(open('boolbits.json'))['boolvars'])
summ=json.load(open('atoms/summary.json')); inputs=set(summ['inputs'])
gate_out=set(t for t,_,_ in gates)
# free inputs = never a gate target
freeinputs=set(v for v in range(NVARS) if v not in gate_out)
# handles = free inputs that are NOT boolean bits (data + absorbers)
handles=sorted(freeinputs - boolbits)
hset=set(handles)
# activation: set ACT bit=1
val0[ACT]=1
for b in freeinputs: pinned[b]=True
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
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[t]][1]),'<r>','eval') for t in order]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
ecode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
# which handles appear in which equations (to build Jacobian only over relevant handles)
eq_handles=[]
for L in lines:
    vs=set(int(m) for m in VAR.findall(L))
    eq_handles.append(vs & hset)
ns={'__builtins__':{}}
val=val0[:]
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[t if False else k],ns)
def resid_all():
    ns['v']=val
    return [eval(c,ns) for c in ecode]
# SNF integer solver (from raw_solve)
def snf_solve(Am,b):
    m=len(Am); nn=len(Am[0]) if Am else 0
    if m==0: return [0]*nn
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
        if bb[i]%d!=0: return None
        y[i]=bb[i]//d
    for i in range(m):
        if all(z==0 for z in M[i]) and bb[i]!=0: return None
    return [sum(V[i][j]*y[j] for j in range(nn)) for i in range(nn)]

forward()
for it in range(12):
    r=resid_all()
    F=[i for i,x in enumerate(r) if x!=0]
    print(f"iter {it}: {len(F)} failing", flush=True)
    if not F: 
        print("SOLVED!"); break
    # relevant handles
    H=sorted(set().union(*[eq_handles[i] for i in F])) if F else []
    hidx={h:j for j,h in enumerate(H)}
    # Jacobian by finite diff
    base=[r[i] for i in F]
    J=[[0]*len(H) for _ in F]
    for j,h in enumerate(H):
        old=val[h]; val[h]=old+1; forward(); rr=resid_all(); val[h]=old
        for ri,i in enumerate(F): J[ri][j]=rr[i]-base[ri]
    forward()
    # solve J*d = -base
    d=snf_solve(J,[-x for x in base])
    if d is None:
        print("  no integer solution to linearized system"); break
    for j,h in enumerate(H): val[h]+=d[j]
    forward()
r=resid_all(); F=[i for i,x in enumerate(r) if x!=0]
print(f"FINAL: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail): {F[:30]}")
if len(F)<26:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('newton_solved.json','w'))
    print("saved newton_solved.json")
