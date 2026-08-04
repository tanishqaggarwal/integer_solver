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
d3=H.loadd('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime3.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d3); H.forward()
# set x_24468, x_18956
x13682=H.val[13682]; x34243=H.val[34243]; x37892=H.val[37892]; x32237=H.val[32237]
d4=dict(d3)
d4[24468]=x13682+12354891*x34243
d4[18956]=x37892+x32237
setfree(d4); H.forward()
F=sorted(H.fails())
print(f"fails: {len(F)}: {F}")
from collections import defaultdict
atom_eqs=defaultdict(list)
for e in F:
    for ai in eq2atoms.get(e,[]):
        if av(ai)!=0: atom_eqs[ai].append(e)
print(f"nonzero atoms: {sorted(atom_eqs.keys())}")
for ai in sorted(atom_eqs):
    print(f"  atom#{ai} n_eq={atoms[ai]['n_eq']} #fails={len(atom_eqs[ai])}: {atoms[ai]['repr'][:70]}")
# M2 detail
L2=H.val[25739]
print(f"\nL2 % p == 0: {L2%p==0}, L2/p mod 6672769 = {(L2//p)%6672769}")
json.dump({f"x_{k}":str(v) for k,v in d4.items()},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime4.json','w'))
