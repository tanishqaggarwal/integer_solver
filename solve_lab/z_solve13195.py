import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
pins=json.load(open('pinrec.json'))
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def mfails():
    return [i for i,c in enumerate(H.eqcode) if eval(c,ns)%p!=0]
print("start x_13195=1:",len(mfails()),"mod-p fails")
# Phase 1: set pins for active bits (targets = CONST) where target is a free input
active_bits=[b for b in set(r[1] for r in pins) if H.val[b]%p==1]
print("active bits:",active_bits)
set_targets=0
for atomidx,sel,target,const,coef,handle in pins:
    if H.val[sel]%p==1:  # active pin -> target should equal const (mod handle)
        if target in H.freeinp:
            setf(target, const%p); set_targets+=1
print(f"set {set_targets} free pin-targets to CONST")
print("after pins:",len(mfails()),"fails")
# Phase 2: zero gadget products. Find atoms A - x*y where target A is 0 -> set free factor to 0
# gadget product vars from diagnosis
for pv in [30175]:  # x_30175=x_26188*x_11075 ; zero via a free ancestor
    anc=[a for a in H.anc.get(pv,set()) if a in H.freeinp]
    # set the free ancestor(s) minimally: try setting each of x_26188,x_11075 factors
# x_30175=x_26188*x_11075: set x_11075=0 if free-reachable
for fac in [26188,11075]:
    fa=[a for a in H.anc.get(fac,{fac}) if a in H.freeinp]
    print(f"x_{fac} free anc: {fa[:5]} (self free={fac in H.freeinp})")
