import heal_harness as H, json
p=H.p
# Load all atoms, find those involving control vars directly
controls={14853,12186,16742,24908,17601,29322,3558}
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        atoms.append(json.loads(line))
print(f"total atoms {len(atoms)}")
for cv in [14853,12186,16742]:
    print(f"\n=== atoms containing x_{cv} directly ===")
    cnt=0
    for ai,a in enumerate(atoms):
        varsin=set()
        for vl,c in a['poly']:
            varsin|=set(vl)
        if cv in varsin:
            cnt+=1
            if cnt<=25:
                print(f"  atom#{ai} neq={a['n_eq']} : {a['repr'][:90]}")
    print(f"  total atoms with x_{cv}: {cnt}")
