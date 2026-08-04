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
d0=H.loadd('best/new_instance_partial_39013.json')
def setb():
    for v in range(H.NVARS): H.val[v]=d0.get(v,0)
setb(); H.forward()
F0=set(H.fails())
# perturb x_22162 by residue-sized amount, check 9648
for delta in [10**18, (10**60)]:
    setb(); H.val[22162]=d0.get(22162,0)+delta
    H.forward()
    F=set(H.fails())
    nb=sorted(F-F0)
    print(f"x_22162 += {delta}: breaks={nb}")
    for e in nb:
        nz=[ai for ai in eq2atoms.get(e,[]) if av(ai)!=0]
        print(f"   eq{e}: nonzero atoms {[(ai,atoms[ai]['repr'][:55]) for ai in nz]}")
# perturb x_30213, check atom29375 and 602
for delta in [10**18, (10**60)]:
    setb(); H.val[30213]=d0.get(30213,0)+delta
    H.forward()
    F=set(H.fails())
    print(f"x_30213 += {delta}: breaks={sorted(F-F0)}, atom29375={av(29375)%p}, atom602={av(602)%p}, atom29371={av(29371)%p}")
