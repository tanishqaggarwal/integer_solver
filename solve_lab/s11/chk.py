import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
import resp as R
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix2_round.json'))
BR=[a for a in range(L.NA) if R.av(v,a)!=0]
C=R.candidates(v,BR,2)
for u in [14853,24548]:
    print(f"x{u} in candidates: {u in C}")
    d1,gb1,n1=R.response(v,u,1); d2,gb2,n2=R.response(v,u,2)
    print(f"   d1 keys={sorted(d1)[:10]} n1={n1}  affine={all(d2.get(a,0)==2*d1.get(a,0) for a in set(d1)|set(d2))}")
    bad=[a for a in set(d1)|set(d2) if d2.get(a,0)!=2*d1.get(a,0)]
    print(f"   non-affine rows: {bad[:10]}")
    for a in bad[:3]:
        print(f"      a{a}: d1={d1.get(a,0)} d2={d2.get(a,0)}")
