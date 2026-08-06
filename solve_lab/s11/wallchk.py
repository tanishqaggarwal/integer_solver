import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
import resp as R
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix7_29539_7930.json'))
BR=[a for a in range(L.NA) if R.av(v,a)!=0]
print("state fix7_29539_7930 broken:",BR)
for u in [30317,5101,21889,32405,5146,32017,2936,26789,25739,37758]:
    print(f"  x{u} {'FREE' if u not in L.definer else 'DEF a%d'%L.definer[u]} val={v[u]} atoms={len(L.var_atoms[u])}")
print()
for u in [30317,5146,2936]:
    if u in L.definer: continue
    d1,_,n=R.response(v,u,1)
    print(f"  x{u}: ripple touches {n} vars; changes checks {sorted(d1)[:12]}")
