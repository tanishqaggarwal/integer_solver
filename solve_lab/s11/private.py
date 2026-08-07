"""The 6 knob atoms of the 15-equation region: do they have PRIVATE handles
   (free variables occurring in exactly one atom)?  If so the equation-space fix is confined."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
KNOBS = [26719, 26720, 26721, 26722, 26723, 28437]

v = [0] * L.NVARS
for k, x in json.load(open(os.path.join(HERE, 'data', 'finish3.json'))).items():
    v[int(k)] = int(x)
fw.forward(v)
av = L.all_atom_values(v)
S = sorted(L.failing_eqs(av))
print(f"failing equations: {len(S)}")


def fmt(a, lim=120):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


steps = {}
for a in KNOBS:
    print(f"a{a} [{len(L.atom2eq.get(a,{}))} eqs] val={'0' if av[a]==0 else 'nz'}: {fmt(a)}")
    try:
        hs, base = deep.handles(v, a, locked=set())
    except Exception:
        hs = []
    priv = [(t, d) for t, d in hs if len(L.var_atoms[t]) == 1]
    print(f"    handles: {[(t, len(L.var_atoms[t])) for t, _ in hs][:8]}   PRIVATE: {[t for t,_ in priv]}")
    if priv:
        steps[a] = priv[0]

print()
print("knobs with a private handle:", sorted(steps))
if len(steps) >= 1:
    ks = sorted(steps)
    # solve  sum_e co[e][a]*(av[a] + delta_a * k_a)  = -(fixed part)  for all e in S
    M = []
    rhs = []
    for e in S:
        mult, sq, co = L.eq_atoms[e]
        row = [co.get(a, 0) * steps[a][1] for a in ks]
        fixed = sum(c * av[a] for a, c in co.items())
        M.append(row)
        rhs.append(-fixed)
    x = solve_int(M, rhs)
    print("confined integer solution over PRIVATE handles:", "FOUND" if x else "NONE")
    if x:
        for j, a in enumerate(ks):
            t = steps[a][0]
            v[t] += x[j]
        fw.forward(v)
        f = L.failing_eqs(L.all_atom_values(v))
        b = fw.bad_checks(v)
        print(f"AFTER: failing={len(f)} score={L.NEQ-len(f)} bad_checks={len(b)} {b[:10]}")
        if len(f) < 15:
            import sys as _s
            _s.set_int_max_str_digits(200000)
            json.dump({('x_%d' % i): str(v[i]) for i in range(L.NVARS)},
                      open(os.path.join(HERE, 'data', 'private_named.json'), 'w'))
            print("saved data/private_named.json")
