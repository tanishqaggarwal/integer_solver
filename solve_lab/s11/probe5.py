import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
def fmt(a, lim=600):
    parts=[]
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]),kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not m) else ('+' if c==1 else '-')
        if m: s+='*'.join('x%d'%u for u in m)
        parts.append(s)
    return ' '.join(parts)[:lim]
for a in [36185, 19297, 19299, 30984, 40812, 16632, 16634, 30976, 30978]:
    print(f"=== a{a} eqs={len(L.atom2eq.get(a,{}))} out={L.atom_out.get(a)}")
    print("   ", fmt(a))
