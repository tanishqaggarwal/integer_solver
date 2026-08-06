import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix2_round.json'))
def sh(x):
    s=str(x); return s if len(s)<22 else s[:9]+'..'+s[-5:]+f'<{len(s)}d>'
for a in [40826,41512,7930,29539]:
    Pp=L.polys[a]
    print(f"=== a{a}  out={L.atom_out.get(a)}  nterms={len(Pp)}  val={sh(atomval(a,v))}  val%p={sh(atomval(a,v)%P)}")
    print("    eqs:", sorted(L.atom2eq.get(a,{})))
    if len(Pp)<40:
        for m,c in Pp.items():
            print(f"      {c} * {'*'.join('x%d'%u for u in m)}   vals={[sh(v[u]) for u in m]}")
    else:
        print(f"      (large: {len(Pp)} monomials, vars={len(L.avars[a])})")
    print()
