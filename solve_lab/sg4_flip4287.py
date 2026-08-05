import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V0=H.val[:]
baseF=set(H.fails())
print("baseline fails:",len(baseF), sorted(baseF))
# record key vars
key=[20492,37158,16738,36065,10878,22542,36760,9062,21279,19892,25297,19964,2099,4432,7068]
rec0={t:V0[t] for t in key}
# STEP 1: flip x_4287 only
H.val[4287]=1
H.forward()
V1=H.val[:]
print("\n=== after x_4287=1 only ===")
for t in key:
    if V1[t]!=rec0[t]:
        d=V1[t]-rec0[t]
        print(f"x_{t}: {rec0[t] if abs(rec0[t])<1e12 else 'BIG'} -> changed by {d if abs(d)<1e12 else str(d)[:10]+'..'} (mod p {d%p})")
F1=set(H.fails())
print("fails after flip:",len(F1), " newly broke:",sorted(F1-baseF)[:25])
# which free inputs are in the newly broken eqs
newbroke=F1-baseF
fv=set()
for i in newbroke: fv|= (H.eqvars[i]&H.freeinp)
print("free inputs in newly-broken eqs:",sorted(fv)[:40])
