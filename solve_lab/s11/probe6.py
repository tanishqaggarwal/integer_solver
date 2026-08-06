import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
def fmt(a, lim=330):
    parts=[]
    for mm,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]),kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not mm) else ('+' if c==1 else '-')
        if mm: s+='*'.join('x%d'%u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]
for a in [26719,26721,26723,26733,28438,32342,36185]:
    print(f"=== a{a} eqs={len(L.atom2eq.get(a,{}))}: {fmt(a)}")
    print("    free:", [(u,len(L.var_atoms[u])) for u in sorted(L.avars[a]) if L.definer.get(u) is None])
