import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
CTRL = [5096, 21589, 14515, 19750, 33708, 31339, 16441, 22917, 13222, 14681, 28486, 38667]


def fmt(a, lim=190):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


for c in CTRL:
    print(f"=== x{c}: in {len(L.var_atoms[c])} atoms")
    for a in L.var_atoms[c]:
        out = L.atom_out.get(a)
        tag = ('GATE->' + str(out)) if out is not None else 'CHECK'
        print(f"    a{a} [{tag}, {len(L.atom2eq.get(a,{}))} eqs]: {fmt(a)}")
