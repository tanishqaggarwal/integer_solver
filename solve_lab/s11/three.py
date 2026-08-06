import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, uv01, uv01full
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def fmt(a, lim=260):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


rng = random.Random(11)
BITS = (490, 91)
v = uv01.state(BITS)
uv01full.structural(v)
LOCK = set(BITS) | {19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742,
                    22162, 30213, 8386, 21868} | set(uv01full.MIRROR_CTRL) | {t for _, t in uv01full.LINK}
for a, t in uv01full.LINK:
    x = fw.solve_lin(a, t, v)
    if x is not None:
        v[t] = x
        fw.forward(v)
uv01full.mirror(v, rng)
uv01full.structural(v)
for a in fw.bad_checks(v):
    if a not in [x for x, _ in uv01full.LINK]:
        uv01full.close_one(v, a, LOCK)
bad = fw.bad_checks(v)
f = L.failing_eqs(L.all_atom_values(v))
print(f"state: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}  {bad}")
print(f"mirror: x3719/P int? {v[3719]%P==0} x25118/P int? {v[25118]%P==0}")
for a in bad:
    val = fw.evalpoly(L.polys[a], v)
    fr = [(u, NAT[u]) for u in sorted(L.avars[a]) if L.definer.get(u) is None]
    print(f"\na{a} [{len(L.atom2eq.get(a,{}))} eqs]  val%P==0: {val%P==0}")
    print(f"   {fmt(a)}")
    print(f"   free: {fr}")
    try:
        hs, base = deep.handles(v, a, locked=LOCK)
        print(f"   cone handles: {[(t, NAT[t], base % d == 0) for t, d in hs][:10]}")
    except Exception as e:
        print("   handles err", e)
json.dump({str(i): v[i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'three.json'), 'w'))
