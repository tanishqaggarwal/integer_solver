"""Regional (equation-space) analysis of the current defect placement.
   S = failing equations.  Knobs = atoms whose ENTIRE equation footprint lies inside S
   (moving them changes nothing outside).  How many of |S| can be recovered?"""
import sys, os, json, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))


def analyse(name):
    v = [0] * L.NVARS
    for k, x in json.load(open(os.path.join(HERE, 'data', name))).items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
    fw.forward(v)
    av = L.all_atom_values(v)
    S = L.failing_eqs(av)
    nz = [a for a, x in enumerate(av) if x != 0]
    print(f"=== {name}: failing={len(S)} nonzero atoms={len(nz)} {nz}")
    # atoms living in these equations
    region = set()
    for e in S:
        region |= set(L.eq_atoms[e][2])
    knobs = [a for a in region if set(L.atom2eq.get(a, {})).issubset(set(S))]
    print(f"    region atoms={len(region)}  knobs (footprint inside S)={len(knobs)}: {sorted(knobs)[:20]}")
    # build M over S x knobs ; rhs = -(contribution of the pinned atoms)
    Sl = sorted(S)
    M = []
    rhs = []
    for e in Sl:
        mult, sq, co = L.eq_atoms[e]
        if sq:
            # equation is mult*(sum)^2 ; vanishes iff sum == 0
            pass
        row = [co.get(a, 0) for a in knobs]
        fixed = sum(c * av[a] for a, c in co.items() if a not in set(knobs))
        M.append(row)
        rhs.append(-fixed)
    x = solve_int(M, rhs) if knobs else None
    print(f"    exact integer solve over knobs: {'FOUND -> all %d recoverable' % len(Sl) if x else 'none'}")
    if not x and knobs:
        # largest recoverable subset (greedy by dropping rows)
        best = 0
        for drop in range(0, min(len(Sl), 9)):
            found = False
            for combo in itertools.combinations(range(len(Sl)), drop):
                keep = [i for i in range(len(Sl)) if i not in combo]
                if solve_int([M[i] for i in keep], [rhs[i] for i in keep]) is not None:
                    best = len(keep)
                    found = True
                    break
            if found:
                print(f"    max recoverable = {best} of {len(Sl)} (drop {drop}) -> "
                      f"failing would be {len(Sl)-best}, score {L.NEQ-(len(Sl)-best)}")
                break
    return len(S)


for nm in ['finish3.json', 'closehit2.json', 'three.json']:
    try:
        analyse(nm)
    except Exception as e:
        print(nm, 'error', e)
