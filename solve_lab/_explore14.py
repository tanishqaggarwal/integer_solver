import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
d=H.loadd('best/new_instance_partial_39013.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
# baseline
setfree(d); H.forward()
base=[H.val[v] for v in range(H.NVARS)]
# regime1
x24908=H.val[24908]
d2=dict(d); d2[14853]=d[12186]; d2[16742]=x24908
setfree(d2); H.forward()
reg=[H.val[v] for v in range(H.NVARS)]
# term-by-term of atom42851
for ai in [42851,44270]:
    print(f"\n===== atom#{ai} term-by-term (baseline -> regime1) =====")
    poly=atoms[ai]['poly']
    for vl,c in poly:
        tb=c
        for v in vl: tb*=base[v]
        tr=c
        for v in vl: tr*=reg[v]
        if tb%p!=tr%p:
            print(f"  CHANGED term coeff={c} vars={vl}: {tb%p} -> {tr%p}")
    # which vars changed?
    vs=set()
    for vl,c in poly: vs|=set(vl)
    changedvars=[v for v in sorted(vs) if base[v]!=reg[v]]
    print(f"  vars in atom that changed value: {changedvars}")
