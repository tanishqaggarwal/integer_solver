import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']; eq2atoms=C['eq2atoms']
def atomval(ai):
    a=atoms[ai]; s=0
    for vl,c in a['poly']:
        t=c
        for v in vl: t*=H.val[v]
        s+=t
    return s
d=H.loadd('best/new_instance_partial_39013.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d); H.forward()
x24908=H.val[24908]
d2=dict(d); d2[14853]=d[12186]; d2[16742]=x24908
setfree(d2); H.forward()
F=sorted(H.fails())
print(f"regime1: {len(F)} fail")
# Map each fail to its nonzero atoms; collect the union with multiplicity of eqs
from collections import defaultdict
atom_eqs=defaultdict(list)
for e in F:
    for ai in eq2atoms.get(e,[]):
        if atomval(ai)!=0:
            atom_eqs[ai].append(e)
print(f"\nDistinct nonzero atoms across all {len(F)} fails: {len(atom_eqs)}")
for ai in sorted(atom_eqs):
    a=atoms[ai]
    print(f"  atom#{ai} n_eq={a['n_eq']} appears_in_fails={atom_eqs[ai]}\n      {a['repr'][:110]}")
# which fails have NO nonzero atom explained above? (should be none)
