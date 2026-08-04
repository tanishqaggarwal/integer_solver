import heal_harness as H
import json
p=H.p
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
# all atoms containing x_7068 or x_4432
atoms=[]; reprs=[]; aeqs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line)
        atoms.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr','')); aeqs.append(dd.get('eqs',[]))
def atomvars(i):
    vs=set()
    for m,c in atoms[i]: vs.update(m)
    return vs
for target in [7068,4432]:
    print(f"=== atoms containing x_{target} (value={H.val[target]%p}) ===")
    for i in range(len(atoms)):
        if target in atomvars(i):
            print(f"  atom {i}: {reprs[i]}  | eqs={aeqs[i]}")
# also equations containing them
print("\n=== equations containing x_7068 ===", [i for i,vs in enumerate(H.eqvars) if 7068 in vs])
print("=== equations containing x_4432 ===", [i for i,vs in enumerate(H.eqvars) if 4432 in vs])
# what are x_2099 and x_19964 (targets)
print("\nx_2099 =", H.val[2099])
print("x_19964 =", H.val[19964])
print("x_7068 =", H.val[7068])
print("x_4432 =", H.val[4432])
print("G1 = x_7068 - x_2099 =", H.val[7068]-H.val[2099])
print("G2 = x_4432 - x_19964 =", H.val[4432]-H.val[19964])
print("G1 mod p =", (H.val[7068]-H.val[2099])%p)
print("G2 mod p =", (H.val[4432]-H.val[19964])%p)
