import os,sys,json,itertools
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
# recompute feasible set quickly (reuse logic, but just re-derive forced bits)
# From prior run: feasible all have 2081=0,4287=0,5910=0,11368=0,13195=1,17406=1? verify count 432
bits10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
# Rather than recompute, let me measure GLOBAL damage of candidate pattern
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F0=set(H.fails()); print("agentA fails:",len(F0))
# candidate: 2081=0,13195=1,17406=1, keep 4287=0, rest of 10 =0 (agentA already 0 for these), other 246 unchanged
cand={2081:0,4287:0,5910:0,11368:0,13195:1,17406:1,18022:0,22562:0,23751:0,28005:0}
for b,v in cand.items(): H.val[b]=v
H.forward()
F1=set(H.fails())
print(f"candidate pattern (2081=0,13195=1,17406=1) total fails: {len(F1)}")
# which loads now nonzero?
for nm,v in [('L1',11150),('L2',25739),('L3',37758),('x_15298',15298),('x_2081bit',2081)]:
    print(f"  {nm}=x_{v}: {H.val[v]%p if v!=15298 else H.val[v]}")
# Now solve local 27 eqs and see global damage
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
LOCAL=FAILS11+RIPPLE16
Eanc=set()
for i in LOCAL:
    for v in H.eqvars[i]: Eanc|=H.anc.get(v,{v})
cont=sorted(v for v in Eanc if v in H.freeinp and v not in bits10)
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
setf(7068,H.val[2099]%p); setf(4432,H.val[19964]%p)
def resid(): return [eval(H.eqcode[i],ns)%p for i in LOCAL]
b0=resid(); m=len(LOCAL); n=len(cont)
Jt=[]
bc={v:H.val[v] for v in cont}
for v in cont:
    o=bc[v]; setf(v,o+1); r1=resid(); setf(v,o)
    Jt.append([(r1[k]-b0[k])%p for k in range(m)])
# solve J * delta = -b0 (least-norm) via flint: augment and solve
MJ=flint.fmpz_mod_mat(m,n,ctx)
for k in range(m):
    for j in range(n): MJ[k,j]=Jt[j][k]
bb=flint.fmpz_mod_mat(m,1,ctx)
for k in range(m): bb[k,0]=(-b0[k])%p
try:
    sol=MJ.solve(bb)
    solvable=True
except Exception as e:
    solvable=False; print("solve failed:",e)
print("local linear system solvable:",solvable)
if solvable:
    for j,v in enumerate(cont): setf(v, (bc[v]+int(sol[j,0]))%p)
    # global fails after applying local solution mod p
    H.forward()  # full re-forward
    F2=set(H.fails())
    print(f"after applying local mod-p solve: total fails={len(F2)}")
    print(f"  local 27 still failing: {len(set(LOCAL)&F2)}")
    print(f"  new global fails: {len(F2-set(LOCAL))}")
