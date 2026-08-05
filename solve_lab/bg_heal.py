#!/usr/bin/env python3
"""Background heal in x_9062=1 quadrant: apply gadget decomposition + iterative try_set heal."""
import heal_harness as H
import json,ast,re,time
from collections import defaultdict
p=H.p
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
d=json.load(open('sy_regime11_39018.json'))
for v in H.freeinp: H.val[v]=0
for k,vv in d.items(): H.val[int(k[2:])]=int(vv)
# apply gadget decomposition: x_14865=x_12553, x_31861=x_6418, x_33168=0
H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0
H.forward()
best=len(set(H.fails()))
print(f"start (gadget-fixed): {best} fails",flush=True)
# descendants per free (for incremental)
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
# iterative: for each failing eq, find a free knob (via its vars) that reduces total fails (greedy, exact)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
def fails(): return set(i for i,c in enumerate(H.eqcode) if eval(c,ns)!=0)
import random
random.seed(1)
F=fails()
t0=time.time()
# candidate knobs: free inputs in the failing eqs' cones + private frees
pool=set()
for i in F:
    for w in H.eqvars[i]:
        if w in H.freeinp: pool.add(w)
        pool|=(H.anc.get(w,set())&H.freeinp)
pool=list(pool)
print(f"pool: {len(pool)} free knobs",flush=True)
# greedy: try setting each free knob to a value that zeros one failing eq (via exact solve of that eq in that var)
for it in range(2000):
    if time.time()-t0>500: print("time budget",flush=True); break
    F=fails()
    if len(F)<best:
        best=len(F); json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)},open('bg_best.json','w'))
        print(f"  it{it}: NEW BEST {39033-best}/39033",flush=True)
        if best==0: print("*** SOLVED ***",flush=True); json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)},open('BG_SOLVED.json','w')); break
    # pick a failing eq, a var in it, solve for that var (linear approx)
    i=random.choice(list(F))
    cand=[w for w in H.eqvars[i] if w in H.freeinp]
    if not cand: continue
    w=random.choice(cand)
    base=eval(H.eqcode[i],ns)
    H.val[w]+=1
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    slope=eval(H.eqcode[i],ns)-base
    if slope!=0 and base%slope==0:
        H.val[w]+= -base//slope -1
    else:
        H.val[w]-=1
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    nf=len(fails())
    if nf>best+3:  # revert if much worse
        pass
print(f"final best: {39033-best}/39033",flush=True)
