"""Cheap gate atoms: breaking one frees its output variable as a NEW control, at the cost of
   the equations it occupies.  Which cheap ones give control over the mirror (x_3719, x_25118)?"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, uv01
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))


def fmt(a, lim=110):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


gp = sorted(((len(L.atom2eq.get(a, {})), a) for a in range(L.NA) if L.atom_out.get(a) is not None))
CHEAP = [(n, a) for n, a in gp if n <= 8]
print(f"gate atoms in <=8 equations: {len(CHEAP)}")
for n, a in CHEAP[:20]:
    out = L.atom_out.get(a)
    ov = out[1] if isinstance(out, tuple) else out
    print(f"  a{a} [{n} eqs] -> x{ov}: {fmt(a)}")

# baseline state in the best branch
v0 = uv01.state((490, 91))
FW_ORDER = fw.ORDER


def forward_skip(v, skip):
    for comp in FW_ORDER:
        if len(comp) == 1:
            u = comp[0]
            if L.definer[u] in skip:
                continue
            x = fw.solve_lin(L.definer[u], u, v)
            if x is not None:
                v[u] = x
        else:
            for _ in range(40):
                ch = False
                for u in comp:
                    if L.definer[u] in skip:
                        continue
                    x = fw.solve_lin(L.definer[u], u, v)
                    if x is not None and x != v[u]:
                        v[u] = x
                        ch = True
                if not ch:
                    break
    return v


print("\nwhich cheap gates give MOD-P control of the mirror when broken?")
hits = []
for n, a in CHEAP:
    out = L.atom_out.get(a)
    ov = out[1] if isinstance(out, tuple) else out
    v = v0[:]
    forward_skip(v, {a})
    b3, b25 = v[3719] % P, v[25118] % P
    v[ov] = v[ov] + 1
    forward_skip(v, {a})
    d3 = (v[3719] % P - b3) % P
    d25 = (v[25118] % P - b25) % P
    if d3 or d25:
        hits.append((n, a, ov, d3 != 0, d25 != 0))
        print(f"  a{a} [{n} eqs] -> x{ov}: moves x3719={d3!=0} x25118={d25!=0}", flush=True)
print(f"\n{len(hits)} cheap gates give mirror control")
json.dump([[n, a, ov] for n, a, ov, _, _ in hits], open(os.path.join(HERE, 'data', 'cheapgates.json'), 'w'))
