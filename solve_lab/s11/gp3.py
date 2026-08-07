import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<24 else s[:10]+'..'+s[-6:]+f'<{len(s)}d>'
def show(u):
    print(f"=== x{u}  {'FREE' if u not in L.definer else 'DEF by a%d'%L.definer[u]}  val={sh(v[u])}  val%p={'0' if v[u]%P==0 else sh(v[u]%P)}")
    for a in L.var_atoms[u]:
        terms=' + '.join(f"{c}*{'*'.join('x%d'%t for t in m)}" for m,c in L.polys[a].items())
        eqs=sorted(L.atom2eq.get(a,{}))
        print(f"   a{a} out={str(L.atom_out.get(a)):13s} val={sh(atomval(a,v))} eqs={eqs[:6]}{'..' if len(eqs)>6 else ''}")
        print(f"      {terms}")
for u in [9118, 8731, 28730, 9413, 1329, 10903]:
    show(u); print()
