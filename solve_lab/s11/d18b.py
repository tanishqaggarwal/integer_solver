import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','finish3_named.json'))
def sh(x):
    s=str(x); return s if len(s)<18 else s[:7]+'..'+s[-4:]+f'<{len(s)}d>'
for u in [3896,12000,24326,24175,35019,12926,21693,4615,16732,21364,3358,13992,35872]:
    d=L.definer.get(u)
    t=(' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in L.polys[d].items())[:80]) if d is not None else 'FREE'
    print(f"x{u:6d}={sh(v[u]):24s} {'==p' if v[u]==P else ''}  atoms={len(L.var_atoms[u])}  {t}")
print()
M=v[3896]
for nm,val in [('x3896*x12000',v[3896]*v[12000]),('x3896*x12926',v[3896]*v[12926]),
               ('2648967*x3896*x21364',2648967*v[3896]*v[21364])]:
    print(f"{nm}: %p={'0' if val%P==0 else 'nz'}  %8640431={val%8640431}")
print("x3896 %p =", 'ZERO' if v[3896]%P==0 else 'nonzero', "  x3896 % 8640431 =", v[3896]%8640431)
