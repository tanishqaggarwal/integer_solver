import heal_harness as H, json, time
p=H.p
d0=H.loadd('best/new_instance_partial_39013.json')
def setb():
    for v in range(H.NVARS): H.val[v]=d0.get(v,0)
setb()
t0=time.time(); H.forward(); print(f"forward time {time.time()-t0:.3f}s")
F0=set(H.fails())
S0=H.val[35389]%p; T0=H.val[6671]%p
# scan cones of x_24908, x_19083 (the b side) and x_1308,x_23927 (a side, not in S,T cone)
stcone=H.anc[35389]|H.anc[6671]
scan=sorted((H.anc[24908]|H.anc[19083]|H.anc[1308]|H.anc[23927]))
print(f"scanning {len(scan)} frees (b-side + a-side cones)")
delta=1234567890123456789
results=[]
for f in scan:
    setb(); H.val[f]=d0.get(f,0)+delta
    H.forward()
    dS=(H.val[35389]%p-S0)%p; dT=(H.val[6671]%p-T0)%p
    F=set(H.fails()); nb=len(F-F0)
    if dS or dT:
        results.append((f,dS!=0,dT!=0,nb,sorted(F-F0)[:5]))
# report cleanest
results.sort(key=lambda r:r[3])
print("frees moving S or T, sorted by #breaks:")
for f,ds,dt,nb,br in results[:25]:
    print(f"  x_{f}: dS={'Y' if ds else '.'} dT={'Y' if dt else '.'} breaks={nb} {br if nb else ''} {'<<CLEAN' if nb==0 else ''}")
print(f"\ntotal moving S/T: {len(results)}, clean: {sum(1 for r in results if r[3]==0)}")
cleanS=[r[0] for r in results if r[3]==0 and r[1]]
cleanT=[r[0] for r in results if r[3]==0 and r[2]]
print(f"clean S-knobs: {cleanS[:10]}")
print(f"clean T-knobs: {cleanT[:10]}")
