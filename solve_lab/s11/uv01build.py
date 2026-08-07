import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, uv01
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def fmt(a, lim=200):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


def drive(v, ctrl, get, target=0, tries=4):
    """set free var ctrl so that get(v) == target mod P (assumes linear)"""
    for _ in range(tries):
        r = (get(v) - target) % P
        if r == 0:
            return True
        old = v[ctrl]
        v[ctrl] = old + 1
        fw.forward(v)
        s = ((get(v) - target) % P - r) % P
        v[ctrl] = old
        fw.forward(v)
        if s == 0:
            return False
        v[ctrl] = (old + (-r) * pow(s, -1, P)) % P
        fw.forward(v)
    return (get(v) - target) % P == 0


BITS = (490, 91)
v = uv01.state(BITS)
print("initial bad:", fw.bad_checks(v))
for a in fw.bad_checks(v):
    fr = [(u, NAT[u]) for u in sorted(L.avars[a]) if L.definer.get(u) is None]
    print(f"  a{a} [{len(L.atom2eq.get(a,{}))} eqs]: {fmt(a)}")
    print(f"     free: {fr}")
print()
# channel invariants (x15574 = 1)
v[22162] = 0
v[30213] = 0
v[8386] = 0
v[21868] = 0
fw.forward(v)
# a688 : x37892 = x24908 ; control x19750
ok1 = drive(v, 19750, lambda vv: vv[37892], G0)
fw.forward(v)
num = C0 + MM * v[37892]
v[7497] = num // P if num % P == 0 else 0
v[22820] = 0
fw.forward(v)
# a1618 : x13682 = x14853 (free)
v[14853] = (-C0B) % P
fw.forward(v)
d = v[13682] + C0B
v[14393] = 0
v[11436] = (d // P) if d % P == 0 else 0
fw.forward(v)
# a29539 : x1308 == x14853 ; control x14515
ok2 = drive(v, 14515, lambda vv: vv[1308] - vv[14853], 0)
# a26731 : x16742 == x19083 (x16742 free & unused in this channel)
v[16742] = v[19083]
fw.forward(v)
print(f"drives: a688={ok1} a29539={ok2}")
print("  a688=%s a1618=%s a40608=%s" % (fw.evalpoly(L.polys[688], v) == 0,
                                        fw.evalpoly(L.polys[1618], v) == 0,
                                        fw.evalpoly(L.polys[40608], v) == 0))
print("  mirror x3719%%p==0:%s x25118%%p==0:%s" % (v[3719] % P == 0, v[25118] % P == 0))
bad = fw.bad_checks(v)
av = L.all_atom_values(v)
f = L.failing_eqs(av)
print(f"  after structural: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}  {bad}")
json.dump({str(i): v[i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'uv01_struct.json'), 'w'))
