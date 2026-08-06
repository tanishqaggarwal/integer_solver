import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, close3
P = L.P
CFGS = {
    'a only          ': (542,),
    'c only          ': (438,),
    'a,b  (U=1,V=0)  ': (542, 47),
    'c,d  (U=0,V=1)  ': (438, 91),
    'a,b,c,d         ': (542, 47, 438, 91),
    'a,c             ': (542, 438),
}
for name, BITS in CFGS.items():
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    fw.forward(v)
    bad = fw.bad_checks(v)
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    print(f"{name}: U={v[7715]} V={v[34554]} x15298={v[15298]} x34606={v[34606]} x5647={v[5647]} "
          f"| bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}")
    print(f"     bad={bad}")
