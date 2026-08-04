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
F=set(H.fails())
newf=sorted(F-set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]))
# newf may include 11854,29437,32916 (core). remove those
core3={11854,29437,32916}
wiring=[e for e in sorted(F) if e not in core3 and e not in set([2071,4573,7123,7469,13660,15299,16622,17726,21382,22093,25480,25539,28653,31061,32894,34517,34892])]
print("wiring fails to analyze:", wiring)
# For each wiring fail, find nonzero atoms
allnz=set()
for e in wiring:
    nzc=[(ai,atomval(ai)) for ai in eq2atoms.get(e,[]) if atomval(ai)!=0]
    print(f"\neq {e}: {len(eq2atoms.get(e,[]))} atoms, {len(nzc)} nonzero")
    for ai,val in nzc:
        allnz.add(ai)
        print(f"   atom#{ai} n_eq={atoms[ai]['n_eq']} val%p={val%p} : {atoms[ai]['repr'][:80]}")
print("\n=== union of nonzero atoms across wiring fails ===", sorted(allnz))
