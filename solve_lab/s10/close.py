"""S10 step 5: try to close  x_7068 == x_2099  from either side, and see what pins it.

Direction A: move x_7068 down to x_2099 (x_7068 is atom 22229's canonical output).
Direction B: move x_2099 up to x_7068 (via its definer x_2099 = x_37158 + x_25297).
Report, for each, exactly which atoms break and in how many equations.
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
CLEAN = os.path.join(HERE, 'clean_state.json')

def report(tag, seeds, block=()):
    v = L.load(CLEAN)
    changed, steps = L.ripple(v, dict(seeds), block=block)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'\n--- {tag}: seeds={ {k:str(x)[:24] for k,x in seeds.items()} }')
    print(f'    changed {len(changed)} vars; nonzero atoms {len(nz)}; '
          f'failing {len(fail)} -> score {L.NEQ-len(fail)}')
    for a in nz:
        print(f'      a{a:<6} neq={len(L.atom2eq.get(a,{})):<3} {L.atom_src[a][:100]}')
    return v, av, fail

v0 = L.load(CLEAN)
print('clean state: x_7068 =', v0[7068])
print('             x_2099 =', v0[2099])
print('             x_37158 =', v0[37158], ' x_25297 =', v0[25297])
print('             x_642 =', v0[642], ' x_21279 =', v0[21279])

# A: pull x_7068 down onto x_2099
report('A  x_7068 := x_2099', {7068: v0[2099]})

# B: push x_2099 up onto x_7068 by moving its input x_37158
report('B  x_37158 += (x_7068-x_2099)', {37158: v0[37158] + (v0[7068] - v0[2099])})
