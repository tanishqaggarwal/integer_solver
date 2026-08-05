import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
pins=json.load(open('pinrec.json'))
selectors=sorted(set(r[1] for r in pins))
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
# close gaps: channel (1,0). set x_7068=x_2099, x_4432=x_19964
setfree(7068, H.val[2099]); setfree(4432, H.val[19964])
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
targetE=FAILS11+RIPPLE16
def resid_vec():
    return tuple(eval(H.eqcode[i],ns)%p for i in targetE)
R0=resid_vec()
# bit support
support=[]
for b in selectors:
    old=H.val[b]; new=1-old if old in (0,1) else 0
    setfree(b,new)
    R1=resid_vec()
    setfree(b,old)
    if R1!=R0:
        nch=sum(1 for a,bb in zip(R0,R1) if a!=bb)
        support.append((b,nch))
print(f"bits whose flip changes the 27-eq residual: {len(support)}")
for b,nch in support:
    print(f"  x_{b}: changes {nch}/27 residuals  (agentA={base[b]})")
# also restore gaps and get the true bit-support ignoring gap-closing (raw 11 fails)
setfree(7068, base[7068]); setfree(4432, base[4432])
R0=tuple(eval(H.eqcode[i],ns)%p for i in FAILS11)
sup11=[]
for b in selectors:
    old=H.val[b]; new=1-old if old in (0,1) else 0
    setfree(b,new); R1=tuple(eval(H.eqcode[i],ns)%p for i in FAILS11); setfree(b,old)
    if R1!=R0: sup11.append(b)
print(f"\nbits affecting the raw 11 fails: {len(sup11)} -> {sup11}")
