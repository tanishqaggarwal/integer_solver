import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
BITS = (542, 47, 438, 91)
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
TH = {int(k): v for k, v in json.load(open('theta_solveA.json')).items()}
v = [0] * L.NVARS
for b in BITS:
    v[b] = 1
for k, x in TH.items():
    v[k] = x
fw.forward(v)
v[14853] = v[12186]
v[16742] = v[24908]
v[30213] = G0
v[22820] = 0
v[7497] = (C0 + MM * G0) // P
v[22162] = -C0B
v[14393] = 0
v[11436] = 0
fw.forward(v)


def fmt(a, lim=200):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


b = fw.bad_checks(v)
print("bad:", len(b), b)
for a in b:
    fr = [(u, len(L.var_atoms[u])) for u in sorted(L.avars[a]) if L.definer.get(u) is None]
    val = fw.evalpoly(L.polys[a], v)
    print(f"a{a} [{len(L.atom2eq.get(a,{}))} eqs] modp={'0' if val%P==0 else 'nz'}: {fmt(a)}")
    print(f"    free: {fr}")
