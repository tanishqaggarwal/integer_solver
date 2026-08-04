import heal_harness as H, json
p=H.p
d0=H.loadd('best/new_instance_partial_39013.json')
def setb():
    for v in range(H.NVARS): H.val[v]=d0.get(v,0)
setb(); H.forward()
F0=set(H.fails())  # 20 core
S0=H.val[35389]%p; T0=H.val[6671]%p
seeds=sorted(H.anc[35389]|H.anc[6671])  # 97 frees feeding S,T
print(f"scanning {len(seeds)} free inputs in S,T cone")
delta=1234567890123456789
clean=[]
for f in seeds:
    setb(); H.val[f]=d0.get(f,0)+delta
    H.forward()
    dS=(H.val[35389]%p - S0)%p
    dT=(H.val[6671]%p - T0)%p
    F=set(H.fails())
    wiring_breaks=sorted((F-F0))   # new fails beyond the core
    nb=len(wiring_breaks)
    if (dS or dT):
        tag="CLEAN" if nb==0 else f"{nb} breaks"
        clean.append((f,dS!=0,dT!=0,nb))
        if nb<=2:
            print(f"  x_{f}: dS={'Y' if dS else '.'} dT={'Y' if dT else '.'}  breaks={nb} {wiring_breaks[:6]}  {'<<< CLEAN' if nb==0 else ''}")
print(f"\nfrees that move S or T: {len(clean)}")
print(f"clean (0 breaks): {[c[0] for c in clean if c[3]==0]}")
print(f"<=1 break: {[(c[0],c[3]) for c in clean if c[3]<=1]}")
