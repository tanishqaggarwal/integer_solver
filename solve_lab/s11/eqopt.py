"""Equation-space optimum for the (490,91) branch:
   three realisable knob atoms  a26719, a26721, a26723  with lattice steps
   (-8640431*p, -p, -p).  Maximise the number of the 15 equations that vanish."""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))

v = [0] * L.NVARS
for k, x in json.load(open(os.path.join(HERE, 'data', 'finish3.json'))).items():
    v[int(k)] = int(x)
fw.forward(v)
av = L.all_atom_values(v)
S = sorted(L.failing_eqs(av))
KN = [(26719, 24175), (26721, 4615), (26723, 13992)]
deltas = []
for a, t in KN:
    old = v[t]
    v[t] = old + 1
    fw.forward(v)
    a2 = L.all_atom_values(v)
    d = a2[a] - av[a]
    v[t] = old
    fw.forward(v)
    deltas.append(d)
    print(f"a{a} <- x{t}: step {d//P if d % P == 0 else d} * p")

# equation e vanishes iff  sum_a co[e][a]*val(a) == 0
rows = []
rhs = []
for e in S:
    mult, sq, co = L.eq_atoms[e]
    row = [co.get(KN[j][0], 0) * deltas[j] for j in range(3)]
    fixed = sum(c * av[a] for a, c in co.items())
    rows.append(row)
    rhs.append(-fixed)

print(f"\n{len(S)} equations, 3 lattice knobs")
best = None
t0 = time.time()
for drop in range(0, 12):
    found = None
    for combo in itertools.combinations(range(len(S)), drop):
        keep = [i for i in range(len(S)) if i not in combo]
        x = solve_int([rows[i] for i in keep], [rhs[i] for i in keep])
        if x is not None:
            found = (combo, x)
            break
    if found:
        combo, x = found
        print(f"  MAX SATISFIABLE = {len(S)-drop} of {len(S)}  (drop {drop}: eqs {[S[i] for i in combo]})")
        print(f"  k = {[str(t)[:24] for t in x]}   ({time.time()-t0:.0f}s)")
        best = (drop, x)
        break
if best:
    drop, x = best
    for j, (a, t) in enumerate(KN):
        v[t] += x[j]
    fw.forward(v)
    f = L.failing_eqs(L.all_atom_values(v))
    b = fw.bad_checks(v)
    print(f"\nAPPLIED: failing={len(f)} score={L.NEQ-len(f)} bad_checks={len(b)}")
    sys.set_int_max_str_digits(300000)
    json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
              open(os.path.join(HERE, 'data', 'eqopt_named.json'), 'w'))
    print("saved data/eqopt_named.json")
