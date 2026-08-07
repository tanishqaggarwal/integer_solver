import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
for nm,f in [('checkpoint',os.path.join(LAB,'best','new_instance_partial_39026.json')),
             ('fix7 cascade end',os.path.join(HERE,'data','fix7_29539_7930.json'))]:
    v=load_raw(f)
    print(f"{nm}:  x15298={v[15298]}  x5101==p:{v[5101]==P} x32017==p:{v[32017]==P} x26789==p:{v[26789]==P}")
    for u in [11150,25739,37758]:
        print(f"    x{u} = {'0' if v[u]==0 else str(v[u])[:12]+'..'}  divisible by p: {v[u]%P==0}")
