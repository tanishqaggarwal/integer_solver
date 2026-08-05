import heal_harness as H
from collections import defaultdict
import random
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()

# Build incremental-forward: desc_of[free] = list of gate-order indices k whose target depends on free
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]:
        desc_of[w].append(k)
# eqs touched by each free input
eq_touch=defaultdict(set)
var2eq=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var2eq[v].append(i)
# map free -> equations: eqs whose vars include free, or include a gate downstream of free
# gate targets downstream of free:
import sys
def eqs_of_free(w):
    s=set(var2eq[w])
    for k in desc_of[w]:
        t=H.order[k]
        s.update(var2eq[t])
    return s

# total nnz estimate over ALL free inputs
base_val=H.val[:]
def incr_set(w, newv):
    H.val[w]=newv
    for k in desc_of[w]:
        H.val[H.order[k]]=eval(H.gcode[k],{'v':H.val,'__builtins__':{}})
def restore():
    H.val[:]=base_val

# Test linearity on relevant free inputs vs their touched eqs
allrel_free = set()
break23=[3408,3841,4134,4526,5069,7276,15440,15724,15927,21600,22139,22825,27289,27999,28718,29305,31134,31269,32463,33195,36387,36390,38888]
F=[2071, 4573, 7123, 7469, 11854, 13660, 15299, 16622, 17726, 21382, 22093, 25480, 25539, 28653, 29437, 31061, 32894, 32916, 34517, 34892]
def eq_free_anc(i):
    s=set()
    for v in H.eqvars[i]:
        s |= (H.anc[v] if v in H.anc else ({v} if v in H.freeinp else set()))
    return s
for i in F+break23: allrel_free |= eq_free_anc(i)
print("relevant free:", len(allrel_free))

# linearity test: for random relevant free f and a touched relevant eq i, check E(f+1)-E(f) vs E(f+2)-E(f+1)
ns={'__builtins__':{}}
def evaleq(i):
    ns['v']=H.val
    return eval(H.eqcode[i],ns)%p
random.seed(1)
nlin=0; nnonlin=0; checked=0
sample=random.sample(sorted(allrel_free), min(40,len(allrel_free)))
relset=set(F+break23)
for f in sample:
    touched=eqs_of_free(f)&relset
    if not touched: continue
    i=sorted(touched)[0]
    b0=evaleq(i)
    incr_set(f, base_val[f]+1); b1=evaleq(i)
    incr_set(f, base_val[f]+2); b2=evaleq(i)
    restore()
    d1=(b1-b0)%p; d2=(b2-b1)%p
    checked+=1
    if d1==d2: nlin+=1
    else: nnonlin+=1
print(f"linearity check on relevant f/eq pairs: linear={nlin} nonlinear={nnonlin} of {checked}")
