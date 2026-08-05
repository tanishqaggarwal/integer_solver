import os,sys,json,time,itertools
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
bits10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
# free ancestors of RIPPLE16
Ranc=set()
for i in RIPPLE16:
    for v in H.eqvars[i]: Ranc|=H.anc.get(v,{v})
mov=sorted(v for v in Ranc if v in H.freeinp and v not in bits10 and v not in (7068,4432))
print(f"movable free inputs (gaps held, bits excluded): {len(mov)}")
# local cone for forwarding
need=set()
for i in RIPPLE16: need|=H.eqvars[i]
need|={2099,19964}
cone=set(need); gdef_vids={t:H.gates[H.definer[t]][2] for t in H.order if t in H.definer}
ch=True
while ch:
    ch=False; add=set()
    for t in list(cone):
        if t in gdef_vids:
            for u in gdef_vids[t]:
                if u not in cone: add.add(u)
    if add: cone|=add; ch=True
localorder=[k for k,t in enumerate(H.order) if t in cone]
# restrict desc recompute to cone
desc_of=defaultdict(list)
for k in localorder:
    t=H.order[k]
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def fwd():
    for k in localorder: H.val[H.order[k]]=eval(H.gcode[k],ns)
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def resid(): return [eval(H.eqcode[i],ns)%p for i in RIPPLE16]
m=len(RIPPLE16); n=len(mov)
def test_pattern(pat):
    for b,val in zip(bits10,pat): H.val[b]=val
    fwd()
    setf(7068,H.val[2099]%p); setf(4432,H.val[19964]%p)  # close+hold gaps
    b0=resid()
    # jacobian over mov
    cols=[]
    bc={v:H.val[v] for v in mov}
    Jt=[]
    for v in mov:
        o=bc[v]; setf(v,o+1); r1=resid(); setf(v,o)
        Jt.append([(r1[k]-b0[k])%p for k in range(m)])
    MJ=flint.fmpz_mod_mat(m,n,ctx)
    Mb=flint.fmpz_mod_mat(m,n+1,ctx)
    for k in range(m):
        for j in range(n):
            MJ[k,j]=Jt[j][k]; Mb[k,j]=Jt[j][k]
        Mb[k,n]=(-b0[k])%p
    return MJ.rank(), Mb.rank()
t0=time.time()
feas=[]; infeas_dims=defaultdict(int)
allpat=list(itertools.product([0,1],repeat=10))
for pat in allpat:
    rJ,rJb=test_pattern(pat)
    if rJ==rJb: feas.append(pat)
    else: infeas_dims[rJb-rJ]+=1
print(f"tested 1024 in {time.time()-t0:.0f}s")
print(f"FEASIBLE patterns (16 ripple zeroable, gaps held): {len(feas)}")
for pat in feas[:10]: print("   ",pat,[bits10[i] for i in range(10) if pat[i]])
print(f"infeasible obstruction-dim distribution: {dict(infeas_dims)}")
