import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
import fw
def fmt(a, lim=400):
    parts=[]
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]), kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not m) else ('+' if c==1 else '-')
        if m: s+='*'.join('x%d'%u for u in m)
        parts.append(s)
    return ' '.join(parts)[:lim]
v=[0]*L.NVARS
fw.forward(v)
for a in [688, 1618, 23000, 39067, 40608, 41211]:
    print(f"=== a{a}  eqs={len(L.atom2eq.get(a,{}))}  value={fw.evalpoly(L.polys[a],v)}")
    print("   ", fmt(a))
    for u in sorted(L.avars[a]):
        d=L.definer.get(u)
        print(f"      x{u}={v[u]!s:.60}  {'FREE' if d is None else 'def a%d'%d}  natoms={len(L.var_atoms[u])}")
