import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict,Counter
p=H.p
pins=json.load(open('pinrec.json'))
selectors=sorted(set(r[1] for r in pins))
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
# descendant recompute map
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
# only need to check eqs that depend on v; precompute eq->vars? Use eqvars
# Build var->eqs index
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
def fails_touching(v):
    # eqs that could change when v (and its descendants) change
    aff=set(var_eqs[v])
    for k in desc_of[v]:
        t=H.order[k]
        aff.update(var_eqs[t])
    return aff
F0=set(H.fails())
print(f"baseline fails: {len(F0)}")
t0=time.time()
results=[]
for b in selectors:
    old=base[b]; new=1-old if old in (0,1) else (1 if old==0 else 0)
    aff=fails_touching(b)
    # eval affected eqs before
    setfree(b,new)
    # count fails among ALL eqs quickly: only aff can change; others = F0
    changed_fail=set(i for i in aff if eval(H.eqcode[i],ns)!=0)
    unaff_fail=F0 - aff
    total=len(changed_fail)+len(unaff_fail)
    setfree(b,old)
    results.append((total,b))
results.sort()
print(f"scan done in {time.time()-t0:.0f}s")
print("bits whose flip gives <= 40 fails:")
for total,b in results:
    if total<=40: print(f"  x_{b}: {total} fails  (agentA val={base[b]})")
print(f"\ndistribution of fail-counts after single-bit flip:")
c=Counter(total for total,b in results)
for k in sorted(c)[:15]: print(f"  {k} fails: {c[k]} bits")
