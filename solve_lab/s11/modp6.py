import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
v=load_raw(os.path.join(HERE,'data','modp5_out.json'))
v0=load_raw(os.path.join(HERE,'data','fix7_29539_7930.json'))
def sh(x):
    s=str(x); return s if len(s)<18 else s[:8]+'..'+f'<{len(s)}d>'
for a in [688,1618,40608,19299,21617]:
    Pp=L.polys[a]
    t=' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in Pp.items())
    print(f"a{a} out={L.atom_out.get(a)} val={sh(L.evalpoly(Pp,v))} (was {sh(L.evalpoly(Pp,v0))}) nterms={len(Pp)}")
    print(f"    {t[:230]}")
    print(f"    free vars: {[u for u in sorted(L.avars[a]) if u in FREE]}")
