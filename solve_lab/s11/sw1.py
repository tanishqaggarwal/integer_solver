import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
def show(u):
    print(f"### x{u} {'FREE' if u not in L.definer else 'DEF a%d'%L.definer[u]} val={sh(v[u])} in {len(L.var_atoms[u])} atoms")
    for a in L.var_atoms[u]:
        Pp=L.polys[a]
        t=' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in Pp.items())
        print(f"  a{a} out={str(L.atom_out.get(a)):12s} neq={len(L.atom2eq.get(a,{}))} val={sh(L.evalpoly(Pp,v))}  {t[:220]}")
for u in [7075,4287,9062,20434,2081]: show(u); print()
