import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
def fmt(a):
    parts=[]
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]), kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not m) else ('+' if c==1 else '-')
        if m: s+='*'.join('x%d'%u for u in m)
        parts.append(s)
    return ' '.join(parts)
for v in [int(x) for x in sys.argv[1:]]:
    print(f"=== x{v}  definer={L.definer.get(v)}  atoms={L.var_atoms[v]}")
    for a in L.var_atoms[v]:
        print(f"   a{a} out=x{L.atom_out.get(a)} eqs={len(L.atom2eq.get(a,{}))}: {fmt(a)[:220]}")
