import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p

core = {14853, 12186, 16742, 24908}
# chain sources' computed sides and pin targets
targets = {
 'x_2099 (G1 target)': 2099, 'x_19964 (G2 target)': 19964,
 'x_1308 (a-pin)': 1308, 'x_23927 (b-pin)': 23927,
 'x_19083 (d-pin)': 19083, 'x_17601 (c-pin)': 17601,
}
# free-input ancestors of each; check intersection with core (and with each other's primaries)
chainA = [7068,2964,23238,2498,28246,6083,14853]
chainB = [4432,24548,36462,14623,11080,31339]
primaries = set(chainA+chainB)

for name, x in targets.items():
    a = H.anc.get(x, set())
    print(f"{name}: #anc={len(a)}")
    print(f"   ∩ core {sorted(a & core)}")
    print(f"   ∩ chain-primaries {sorted(a & primaries)}")

# Also: does x_12186/x_16742 feed the chain sources? i.e. are core vars ancestors of x_2099 etc?
# anc gives free-input ancestors. core vars 14853,12186,16742 are themselves free inputs.
print("\n=== Are core free-vars ancestors of the targets? (cycle test) ===")
for name, x in targets.items():
    a = H.anc.get(x, set())
    feeders = sorted(a & core)
    print(f"  {name} fed by core vars: {feeders}")

# key consistency conditions residues at fullcore_fix
print("\n=== residues at fullcore_fix ===")
vf = H.loadd('fullcore_fix.json')
for v in H.freeinp: H.val[v]=vf.get(v,0)
H.forward()
def r(x): return H.val[x]%p
for lbl,(x,y) in {'a=x_1308':(14853,1308),'b=x_23927':(12186,23927),
                  'd=x_19083':(16742,19083),'c=x_17601':(24908,17601)}.items():
    print(f"  {lbl}: x_{x}%p={r(x)==r(y)} (match={r(x)==r(y)})  need for pin")
