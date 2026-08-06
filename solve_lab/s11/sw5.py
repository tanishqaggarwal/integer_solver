import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','sw3_out.json'))
for u in [2239,31731,9106,27177,4306,24425,15963,3177,20473,25490,37944,3349,17925,27019,34310,7181,7010,21671,37530]:
    a=L.definer.get(u)
    t=' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in L.polys[a].items()) if a is not None else 'FREE'
    print(f"x{u:6d} %p={'ZERO' if v[u]%P==0 else 'nonzero'}   {t[:110]}")
print()
print("does x8731 move x27177?  d(x27177)/d(x8731):")
v2=list(v); L.ripple(v2,{8731:v[8731]+1})
for u in [27177,4306,2239,31731,9106]:
    print(f"   x{u}: changed={v2[u]!=v[u]}  delta%p={(v2[u]-v[u])%P!=0}")
print()
print("does x9118 move x27177?")
v3=list(v); L.ripple(v3,{9118:v[9118]+1})
for u in [27177,4306,2239,31731,9106]:
    print(f"   x{u}: changed={v3[u]!=v[u]}  delta%p nonzero={(v3[u]-v[u])%P!=0}")
