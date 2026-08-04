import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
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
print("=== BASELINE 39013 congruences mod p ===")
for label,a,b in [("x_14853 vs x_1308",14853,1308),("x_16742 vs x_19083",16742,19083),
                  ("x_24908 vs x_19083",24908,19083),("x_14853 vs x_12186",14853,12186),
                  ("x_1308 vs x_12186",1308,12186)]:
    print(f"  {label}: {(H.val[a]-H.val[b])%p}")
# baseline verifier atoms (should be 0)
print("baseline atom42851,44270,43834,25170,27902 =", av(42851)%p, av(44270)%p, av(43834)%p, av(25170)%p, av(27902)%p)
# Now check: which vars in atom42851 are downstream of x_14853/x_16742/x_12186?
desc=set()
for t in H.order:
    if H.anc[t] & {14853,16742,12186}: desc.add(t)
print("\natom42851 vars downstream of controls:", sorted(set([1308,3114,11360,11597,14515,14853,19750,29967,30163,32572,36808,36977]) & (desc|{14853,16742,12186})))
print("atom44270 vars downstream of controls:", sorted(set([1308,2873,4494,8408,11360,11597,11768,12182,14515,14853,16015,16787,18340,18534,19450,19750,20454,21608,24509,25783,27040,29967,30163,30355,30454,31330,32572,32683,34575,36808,36977]) & (desc|{14853,16742,12186})))
