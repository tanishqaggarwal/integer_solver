import heal_harness as H
from collections import defaultdict
import random
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
ns={'__builtins__':{}}
def eqresid_all():
    ns['v']=H.val
    return [eval(c,ns)%p for c in H.eqcode]

# baseline residuals mod p
r0=eqresid_all()
base_nz=set(i for i,r in enumerate(r0) if r!=0)
print("baseline nonzero mod p:", len(base_nz), sorted(base_nz))

# Apply random perturbations to ALL free inputs, forward, collect breakable eqs
breakable=set()
freelist=sorted(H.freeinp)
for seed in range(6):
    random.seed(seed)
    for v in freelist: H.val[v]=random.randrange(p)
    H.forward()
    ns['v']=H.val
    for i,c in enumerate(H.eqcode):
        if eval(c,ns)%p!=0: breakable.add(i)
    print(f"seed {seed}: cumulative breakable={len(breakable)}")

# restore
for v in H.freeinp: H.val[v]=base[v]
H.forward()
print("TOTAL breakable (non-constructive) eqs:", len(breakable))
# Are the 20 baseline-failing among breakable?
print("baseline-fail subset of breakable:", base_nz<=breakable)
import pickle
pickle.dump(sorted(breakable), open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/breakable.pkl','wb'))
