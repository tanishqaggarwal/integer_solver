import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
def sh(x):
    s=str(x); return s if len(s)<18 else s[:7]+'..'+s[-4:]+f'<{len(s)}d>'
for name,path in [('checkpoint 39026', os.path.join(LAB,'best','new_instance_partial_39026.json')),
                  ('s11 best 39018',  os.path.join(HERE,'data','finish3_named.json'))]:
    v=load_raw(path)
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    F=L.failing_eqs(AV)
    print(f"=== {name}: score {L.NEQ-len(F)}   x15298={v[15298]} x5647={v[5647]} x34606={v[34606]}")
    for a in [19297,19299,30984]:
        Pp=L.polys[a]
        print(f"   a{a} val={sh(AV[a])} vars: " +
              ' '.join(f"x{u}={sh(v[u])}" for u in sorted(L.avars[a])))
    for u in [4007,11150,25739,29804,35605,37758]:
        d=L.definer.get(u)
        t=(' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in L.polys[d].items())[:90]) if d is not None else 'FREE'
        print(f"   x{u}={sh(v[u])}  {t}")
    print()
