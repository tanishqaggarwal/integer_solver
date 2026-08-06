import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v = load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
ATOMS=[22229,22230,35756,35757,35758,35759,35760,35761,35762]
for a in ATOMS:
    Pp=L.polys[a]
    out=L.atom_out.get(a)
    print(f"--- atom {a}  out={out}  nterms={len(Pp)} val={atomval(a,v)}")
    for m,c in list(Pp.items()):
        print("     ", c, m, [v[u] for u in m])
