import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']; eq2atoms=C['eq2atoms']
def av(ai):
    a=atoms[ai]; s=0
    for vl,c in a['poly']:
        t=c
        for v in vl: t*=H.val[v]
        s+=t
    return s
d=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d.get(v,0)
H.forward()
F=sorted(H.fails())
print(f"BASELINE 39013: {len(F)} fail: {F}")
from collections import defaultdict
atom_eqs=defaultdict(list)
for e in F:
    for ai in eq2atoms.get(e,[]):
        if av(ai)!=0: atom_eqs[ai].append(e)
print(f"nonzero atoms: {len(atom_eqs)}")
for ai in sorted(atom_eqs):
    a=atoms[ai]
    vs=set()
    for vl,c in a['poly']: vs|=set(vl)
    freev=sorted(v for v in vs if v in H.freeinp)
    print(f"  atom#{ai} n_eq={a['n_eq']} fails={atom_eqs[ai]}")
    print(f"     {a['repr'][:120]}")
    print(f"     free vars: {freev}")
