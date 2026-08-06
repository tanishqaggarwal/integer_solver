"""General equation-space optimum: at any state, find the atoms whose whole equation
   footprint lies inside the failing set, get their REALISABLE lattice steps from private
   handles, and maximise how many failing equations can be made to vanish."""
import sys, os, json, itertools, time, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(300000)


def analyse(name, maxdrop=8):
    v = [0] * L.NVARS
    for k, x in json.load(open(os.path.join(HERE, 'data', name))).items():
        v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
    fw.forward(v)
    av = L.all_atom_values(v)
    S = sorted(L.failing_eqs(av))
    if not S:
        print(f"{name}: already 0 failing!")
        return 0, v
    region = set()
    for e in S:
        region |= set(L.eq_atoms[e][2])
    knobs = [a for a in region if set(L.atom2eq.get(a, {})).issubset(set(S))]
    # realisable steps: private free handles (var in exactly one atom)
    steps = []
    for a in knobs:
        try:
            hs, _ = deep.handles(v, a, locked=set())
        except Exception:
            hs = []
        pr = [(t, d) for t, d in hs if len(L.var_atoms[t]) == 1 and d]
        if pr:
            steps.append((a, pr[0][0], pr[0][1]))
    print(f"{name}: failing={len(S)} knobs={len(knobs)} realisable={len(steps)} "
          f"{[(a, t) for a, t, _ in steps]}")
    if not steps:
        return len(S), v
    rows, rhs = [], []
    for e in S:
        mult, sq, co = L.eq_atoms[e]
        rows.append([co.get(a, 0) * d for a, t, d in steps])
        rhs.append(-sum(c * av[a] for a, c in co.items()))
    t0 = time.time()
    for drop in range(0, min(maxdrop, len(S)) + 1):
        for combo in itertools.combinations(range(len(S)), drop):
            keep = [i for i in range(len(S)) if i not in combo]
            x = solve_int([rows[i] for i in keep], [rhs[i] for i in keep])
            if x is not None:
                for j, (a, t, d) in enumerate(steps):
                    v[t] += x[j]
                fw.forward(v)
                f = L.failing_eqs(L.all_atom_values(v))
                print(f"    -> satisfiable {len(S)-drop}/{len(S)}; applied: failing={len(f)} "
                      f"score={L.NEQ-len(f)} ({time.time()-t0:.0f}s)")
                return len(f), v
        if time.time() - t0 > 240:
            break
    print(f"    -> no improvement within drop<={maxdrop} ({time.time()-t0:.0f}s)")
    return len(S), v


best = None
for nm in ['closehit2.json', 'three.json', 'finish3.json', 'quad3_hit.json', 'tri7_best.json',
           'uv01_full.json', 'joint_best.json']:
    if not os.path.exists(os.path.join(HERE, 'data', nm)):
        continue
    try:
        f, v = analyse(nm)
    except Exception as e:
        print(nm, 'error', e)
        continue
    if best is None or f < best[0]:
        best = (f, v, nm)
print(f"\nBEST failing={best[0]} score={L.NEQ-best[0]} from {best[2]}")
json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
          open(os.path.join(HERE, 'data', 'eqopt2_named.json'), 'w'))
