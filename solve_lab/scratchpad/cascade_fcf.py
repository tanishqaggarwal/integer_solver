import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

vf = H.loadd('fullcore_fix.json')
for v in H.freeinp: H.val[v] = vf.get(v, 0)
H.forward()
print(f"start fullcore_fix: {len(H.fails())} fails")

def zero_atom(ai):
    poly = A.ATOMS[ai]
    val = A.eval_atom(ai, H.val)
    for vl, c in poly:
        if len(vl) == 1 and vl[0] in H.freeinp:
            f = vl[0]
            if val % c == 0:
                H.val[f] -= val // c
                return f
    return None

seen = {}
for step in range(40):
    sn = [ai for ai in range(A.NATOM) if len(A.ATOM_VARS[ai]) <= 6 and A.eval_atom(ai, H.val) != 0]
    if not sn:
        print(f"step {step}: all small gaps zero")
        break
    for ai in sn: seen.setdefault(ai, A.ATOM_REPR[ai])
    fixed = False
    for ai in sn:
        if zero_atom(ai) is not None: fixed = True
    H.forward()
    if not fixed:
        print(f"step {step}: STUCK. small gaps not linearly fixable: {sn}")
        for ai in sn:
            print(f"    atom {ai}: {A.ATOM_REPR[ai][:80]}")
        break

# terminal residual: ALL nonzero atoms now
nz = A.nonzero_atoms(H.val)
print(f"\nTERMINAL: {len(H.fails())} fails, {len(nz)} nonzero atoms")
freev = set()
for ai, val in nz:
    fv = [v for v in A.ATOM_VARS[ai] if v in H.freeinp]
    freev |= set(fv)
    print(f"  atom {ai}: {A.ATOM_REPR[ai][:70]}  #vars={len(A.ATOM_VARS[ai])} freevars_in_atom={sorted(fv)[:8]}")
print(f"\ntotal distinct gaps healed along the way: {len(seen)}")
print(f"free vars directly in terminal nonzero atoms: {len(freev)}: {sorted(freev)}")
