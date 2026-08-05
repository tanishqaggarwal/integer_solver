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

# classify atoms: "primitive gap" = nonzero atom with few vars (<=6) that is NOT a pure gate def
# a gate-def atom is one whose target is a gate output and appears as x_t - rhs
# We'll just track small nonzero atoms.
def small_nonzero():
    out = []
    for ai in range(A.NATOM):
        if len(A.ATOM_VARS[ai]) <= 6:
            val = A.eval_atom(ai, H.val)
            if val != 0:
                out.append((ai, val))
    return out

# For a small gap atom, identify: which of its vars is FREE (the one to move),
# and set that free var so the atom -> 0 by making free = (rest).
# atom = sum of terms. If exactly one var is free and appears linearly with coeff +-1, solve.
def atom_terms(ai):
    return A.ATOMS[ai]

seen_gaps = {}   # ai -> (repr, freevars)
history = []
for step in range(60):
    sn = small_nonzero()
    if not sn:
        print(f"step {step}: NO small nonzero atoms remaining!")
        break
    # record
    for ai, val in sn:
        fv = [x for x in A.ATOM_VARS[ai] if x in H.freeinp]
        seen_gaps.setdefault(ai, (A.ATOM_REPR[ai], tuple(sorted(A.ATOM_VARS[ai])), tuple(fv)))
    F = H.fails()
    history.append((step, len(sn), len(F), sorted(ai for ai,_ in sn)))
    # fix each small gap: find a free var with linear coeff, set it to zero the atom
    progress = False
    for ai, val in sn:
        poly = A.ATOMS[ai]
        # find free var appearing in a linear (single-var) term with coeff +/-1
        target = None; coeff=None
        for varlist, c in poly:
            if len(varlist) == 1 and varlist[0] in H.freeinp and abs(c) == 1:
                target = varlist[0]; coeff = c; break
        if target is None:
            continue
        # current atom value with target's contribution removed:
        # atom = c*target + rest  => set target = -rest/c to zero it
        rest = val - coeff * H.val[target]
        H.val[target] = (-rest)//coeff if (-rest) % coeff == 0 else None
        if H.val[target] is None:
            H.val[target] = -rest * coeff  # coeff=+-1 so fine
        progress = True
    H.forward()
    if not progress:
        print(f"step {step}: no linearly-fixable small gap (targets are products?)")
        break

print(f"\n=== cascade history (step, #small_nz, #fails, atoms) ===")
for s in history[:40]:
    print(f"  step {s[0]}: small_nz={s[1]} fails={s[2]} atoms={s[3][:12]}")
print(f"\nTotal distinct small gap atoms seen: {len(seen_gaps)}")
for ai in sorted(seen_gaps):
    r, allv, fv = seen_gaps[ai]
    print(f"  atom {ai}: {r[:80]!r}  freevars={fv}")
