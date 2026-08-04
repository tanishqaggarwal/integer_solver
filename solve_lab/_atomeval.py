import heal_harness as H, json
p=H.p
# Build atom evaluator: each atom poly = list of [varlist, coeff]; value = sum coeff*prod(vars)
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
# map eq index -> list of atom indices (via 'eqs' field)
from collections import defaultdict
eq2atoms=defaultdict(list)
for ai,a in enumerate(atoms):
    for e in a['eqs']:
        eq2atoms[e].append(ai)
def atomval(a):
    s=0
    for vl,c in a['poly']:
        t=c
        for v in vl: t*=H.val[v]
        s+=t
    return s
def atomvars(a):
    S=set()
    for vl,c in a['poly']: S|=set(vl)
    return S
# save for reuse
import pickle
pickle.dump({'atoms':atoms,'eq2atoms':dict(eq2atoms)}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','wb'))
print("atoms cached:",len(atoms))
