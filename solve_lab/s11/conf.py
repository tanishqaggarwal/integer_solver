import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fixpoint_out.json'))
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
for a in [21617,12108,19297,19299,30984,36185,40812,25676,33796,40562,42245]:
    Pp=L.polys[a]
    fv=[u for u in sorted(L.avars[a]) if u in FREE]
    t=' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in Pp.items())
    print(f"a{a} out={str(L.atom_out.get(a)):12s} neq={len(L.atom2eq.get(a,{}))} free vars={fv}")
    print(f"    {t[:200]}")
