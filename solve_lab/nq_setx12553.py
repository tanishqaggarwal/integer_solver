import heal_harness as H
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
print(f"baseline fails: {len(F0)}: {sorted(F0)}")
cands=[95633725943103927163511078726842680621987305533354418579536612548753687164483,
       101660107408702603552292264927445892062532150144667839242556857779371564931523]
orig=V[12553]
for c in cands:
    V[12553]=c
    H.forward()
    F=set(H.fails())
    print(f"\nx_12553={c%p if c<p else 'big'} (residue set):")
    print(f"  fails: {len(F)}  x_27177%p={V[27177]%p}  x_31731={V[31731]}")
    print(f"  x_27177={V[27177]}")
    print(f"  fixed vs baseline: {sorted(F0-F)}")
    print(f"  newly broken: {sorted(F-F0)}")
    V[12553]=orig
    H.forward()
