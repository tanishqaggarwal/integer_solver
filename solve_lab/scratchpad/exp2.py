import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()
x2099 = H.val[2099]; x19964 = H.val[19964]; x642 = H.val[642]; x28730 = H.val[28730]
H.val[7068] = x2099 + 7376877 * x642
H.val[4432] = x19964 + x28730
H.forward()
F1 = H.fails()
print(f"after move: {len(F1)} fails")

# find nonzero atoms
nz = A.nonzero_atoms(H.val)
print(f"{len(nz)} nonzero atoms:")
allvars = set()
for ai, val in nz:
    vp = val % p
    print(f"  atom {ai}: {A.ATOM_REPR[ai][:110]!r}")
    print(f"      %p={'0' if vp==0 else 'NZ'} vars={sorted(A.ATOM_VARS[ai])} ineqs={len(A.ATOM_EQS[ai])}")
    allvars |= A.ATOM_VARS[ai]
print(f"\nunion of vars in nonzero atoms ({len(allvars)}): {sorted(allvars)}")
# which of these are free?
freev = [x for x in sorted(allvars) if x in H.freeinp]
print(f"free among them: {freev}")
