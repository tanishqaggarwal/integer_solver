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
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d); H.forward()
F0=set(H.fails())
allatoms=range(len(atoms))
def nzatoms():
    return set(ai for ai in allatoms if av(ai)!=0)
nz0=nzatoms()
print("baseline nonzero atoms:", sorted(nz0))
# Perturb each control individually
for cv in [14853,12186,16742]:
    setfree(d); H.val[cv]=d.get(cv,0)+7777777
    H.forward()
    nz=nzatoms()
    newatoms=sorted(nz-nz0)
    print(f"\nperturb x_{cv}: new nonzero atoms = {newatoms}")
    for ai in newatoms:
        print(f"    atom#{ai} n_eq={atoms[ai]['n_eq']}: {atoms[ai]['repr'][:75]}")
# Is atom25168 a check that pins x_12186? which var does forward define via it?
# Find definer gate for x_23927, x_25758
for v in [23927,25758,12186]:
    if v in H.definer:
        print(f"\nx_{v} definer gate rhs: {H.gates[H.definer[v]][1]}")
    else:
        print(f"\nx_{v} is FREE (no definer)")
