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

CORE=[7930,29539,35759,35760,40826,41512]
for a in CORE:
    eqs=sorted(L.atom2eq.get(a,{}))
    print(f"a{a}  [{len(eqs)} eqs]  out={L.atom_out.get(a)}")
    print("   ", fmt(a)[:300])
    for u in sorted(L.avars[a]):
        d=L.definer.get(u)
        n=len(L.var_atoms[u])
        print(f"      x{u}: {'FREE' if d is None else 'def a%d'%d}, in {n} atoms, in {len(L.var_eqs[u])} eqs")
    print()
