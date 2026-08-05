import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v] = vA.get(v, 0)
H.forward()
F0 = set(H.fails())
print(f"baseline: {len(F0)} fails")
print(f"G1 = x_7068 - x_2099 - 7376877*x_642 = {H.val[7068]-H.val[2099]-7376877*H.val[642]}")
print(f"G2 = x_4432 - x_19964 - x_28730 = {H.val[4432]-H.val[19964]-H.val[28730]}")

# Move x_2099 to equal x_7068 (heal G1) via x_6418 (slope +1 into x_2099); and x_19964 to x_4432 via x_12553
g1 = H.val[7068] - H.val[2099] - 7376877*H.val[642]   # want x_2099 += g1  => x_6418 += g1
g2 = H.val[4432] - H.val[19964] - H.val[28730]
H.val[6418] += g1
H.val[12553] += g2
H.forward()
# verify targets moved as expected
print(f"\nafter x_6418+=g1, x_12553+=g2:")
print(f"  new G1 = {H.val[7068]-H.val[2099]-7376877*H.val[642]}")
print(f"  new G2 = {H.val[4432]-H.val[19964]-H.val[28730]}")
F1 = set(H.fails())
print(f"  fails: {len(F1)}   newbroken={sorted(F1-F0)}   fixed={sorted(F0-F1)}")

# nonzero atoms now
nz = A.nonzero_atoms(H.val)
print(f"  nonzero atoms: {len(nz)}")
for ai, val in nz[:25]:
    fv = [v for v in A.ATOM_VARS[ai] if v in H.freeinp]
    print(f"    atom {ai}: {A.ATOM_REPR[ai][:65]}  freevars={sorted(fv)[:6]}")
