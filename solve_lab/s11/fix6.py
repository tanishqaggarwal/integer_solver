import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix2_round.json'))
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
for u in [24548,14853,7927,25442,1308]:
    print(f"x{u} {'FREE' if u not in L.definer else 'DEF a%d'%L.definer[u]} val={sh(v[u])} in {len(L.var_atoms[u])} atoms:")
    for a in L.var_atoms[u]:
        deg=max(m.count(u) for m in L.polys[a] if u in m)
        print(f"    a{a} out={str(L.atom_out.get(a)):12s} val={sh(atomval(a,v))} deg_u={deg} neq={len(L.atom2eq.get(a,{}))}")
print()
print("a29539:")
for m,c in L.polys[29539].items():
    print("   ",c,'*','*'.join('x%d'%u for u in m), [sh(v[u]) for u in m])
print("a7929:")
for m,c in L.polys[7929].items():
    print("   ",c,'*','*'.join('x%d'%u for u in m), [sh(v[u]) for u in m])
