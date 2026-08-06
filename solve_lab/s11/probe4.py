import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
def fmt(a, lim=300):
    parts=[]
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]), kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not m) else ('+' if c==1 else '-')
        if m: s+='*'.join('x%d'%u for u in m)
        parts.append(s)
    return ' '.join(parts)[:lim]
v=[0]*L.NVARS; v[542]=1; v[438]=1
fw.forward(v)
for a in [13438,13440,14312,36040,36042,41774,688,1618,40608]:
    print(f"=== a{a} eqs={len(L.atom2eq.get(a,{}))} val={str(fw.evalpoly(L.polys[a],v))[:50]}")
    print("   ", fmt(a))
    for u in sorted(L.avars[a]):
        d=L.definer.get(u)
        print(f"      x{u}={str(v[u])[:40]} {'FREE' if d is None else 'def a%d'%d} natoms={len(L.var_atoms[u])}")
print()
print("x37892 =", str(v[37892])[:60], " x13682 =", str(v[13682])[:60])
print("x30213 free?", L.definer.get(30213), " atoms:", L.var_atoms[30213])
print("x22162 free?", L.definer.get(22162), " atoms:", L.var_atoms[22162])
