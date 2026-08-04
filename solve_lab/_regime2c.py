import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
# dependency of atom29373,29375 vars on knobs
knobs={22162,30213,30317,2936,5146}
for ai in [29373,29375]:
    a=atoms[ai]
    vs=set()
    for vl,c in a['poly']: vs|=set(vl)
    print(f"\natom#{ai}: {a['repr']}")
    for v in sorted(vs):
        anc=H.anc.get(v,{v})
        dep=anc & knobs
        print(f"   x_{v}: {'FREE' if v in H.freeinp else 'gate'}  depends on knobs: {sorted(dep) if dep else 'none'}  #anc={len(anc)}")
