import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
for lab,sel in [("U=1,V=0",{542:1}), ("U=0,V=1",{438:1}), ("U=1,V=1",{542:1,438:1})]:
    v=[0]*L.NVARS
    for k,x in sel.items(): v[k]=x
    fw.forward(v)
    b=fw.bad_checks(v)
    av=L.all_atom_values(v); f=L.failing_eqs(av)
    print(f"{lab}: U={v[7715]} V={v[34554]} | x15298={v[15298]} x5647={v[5647]} x34606={v[34606]}")
    print(f"    bad_checks={len(b)} failing_eqs={len(f)} score={L.NEQ-len(f)}  checks={b}")
