import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v = load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
FAILEQ=[12231,12270,12350,14584,18673,22044,29125]
ATOMS=set()
for e in FAILEQ: ATOMS |= set(L.eq_atoms[e][2])
VARS=set()
for a in ATOMS: VARS |= set(L.avars[a])
print("atoms in failing eqs:",len(ATOMS)," vars:",len(VARS))
print()
for u in sorted(VARS):
    outside_atoms = [a for a in L.var_atoms[u] if a not in ATOMS]
    outside_eqs=set()
    for a in outside_atoms: outside_eqs |= set(L.atom2eq.get(a,{}))
    dfn = L.definer.get(u)
    print(f"x{u:6d} {'DEF by a%-6d'%dfn if dfn is not None else 'FREE        ':14s} "
          f"atoms={len(L.var_atoms[u])} outside_atoms={len(outside_atoms)} outside_eqs={len(outside_eqs)} "
          f"{'ISP' if v[u]==P else ''}")
