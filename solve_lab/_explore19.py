import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
for ai in [601,602,1464,1465,18081,18082,18083,18084]:
    print(f"atom#{ai} n_eq={atoms[ai]['n_eq']}: {atoms[ai]['repr'][:130]}")
# Check x_22162, x_30213 roles: what atoms contain them?
print("\n=== atoms containing x_22162 ===")
cnt=0
for ai,a in enumerate(atoms):
    vs=set()
    for vl,c in a['poly']: vs|=set(vl)
    if 22162 in vs:
        cnt+=1
        if cnt<=20: print(f"  atom#{ai} n_eq={a['n_eq']}: {a['repr'][:70]}")
print(f"  total: {cnt}")
print("\n=== atoms containing x_30213 ===")
cnt=0
for ai,a in enumerate(atoms):
    vs=set()
    for vl,c in a['poly']: vs|=set(vl)
    if 30213 in vs:
        cnt+=1
        if cnt<=20: print(f"  atom#{ai} n_eq={a['n_eq']}: {a['repr'][:70]}")
print(f"  total: {cnt}")
