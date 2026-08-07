"""Breaking a GATE atom frees its output variable -> one extra control, at the cost of the
   equations that atom occupies.  Price the gates that would decouple the two collisions:
     a14445 & a34580 both driven only by x_33129   (x_33129 -> x_15111 -> x_20541 -> x_10170)
     a27139 & a33796 both driven only by x_37088
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P


def fmt(a, lim=130):
    parts = []
    for mm, c in sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0])):
        s = ('%+d' % c) if (c not in (1, -1) or not mm) else ('+' if c == 1 else '-')
        if mm:
            s += '*'.join('x%d' % u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]


def chain(v0, target, depth=6):
    """gates on the path from target upward"""
    out = []
    seen = set()
    st = [(target, 0)]
    while st:
        u, d = st.pop()
        if u in seen or d > depth:
            continue
        seen.add(u)
        g = L.definer.get(u)
        if g is None:
            continue
        out.append((u, g, len(L.atom2eq.get(g, {}))))
        for w in L.avars[g]:
            if w != u:
                st.append((w, d + 1))
    return out


print("GATES on the x_33129 -> x_10170 path (breaking one decouples a14445 from a34580):")
for u, g, n in chain(None, 10170, 4):
    print(f"  x{u} := a{g}  [{n} eqs]: {fmt(g)}")
print()
print("GATES on the x_37088 -> x_6858 path (decouples a27139 from a33796):")
for u, g, n in chain(None, 6858, 4):
    print(f"  x{u} := a{g}  [{n} eqs]: {fmt(g)}")
print()
print("cheapest gate atoms overall (by equation count):")
gp = sorted(((len(L.atom2eq.get(a, {})), a) for a in range(L.NA) if L.atom_out.get(a) is not None))
print("  ", gp[:15])
