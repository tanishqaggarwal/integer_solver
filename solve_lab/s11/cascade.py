import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, uv10
P = L.P
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def fmt(a, lim=210):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


BITS = (542, 47)
v = uv10.state(BITS, {})
bad = fw.bad_checks(v)
print(f"BITS={BITS} bad={len(bad)}")
for a in bad:
    fr = [(u, NAT[u]) for u in sorted(L.avars[a]) if L.definer.get(u) is None]
    print(f"a{a} [{len(L.atom2eq.get(a,{}))} eqs] modp={'0' if fw.evalpoly(L.polys[a],v)%P==0 else 'nz'}: {fmt(a)}")
    print(f"     free: {fr}")
