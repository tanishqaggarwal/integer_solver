import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
# leaves: x_8599 = OR(10566,13886,14188,16734) -> U ; x_7304 = OR(5077,13976,36835,16739) -> V
for lab,sel in [("U=0,V=0",{}), ("U=1,V=0",{10566:1}), ("U=0,V=1",{5077:1}), ("U=1,V=1",{10566:1,5077:1})]:
    v=[0]*L.NVARS
    for k,x in sel.items(): v[k]=x
    fw.forward(v)
    b=fw.bad_checks(v)
    av=L.all_atom_values(v); f=L.failing_eqs(av)
    print(f"{lab}: U={v[7715]} V={v[34554]} x15298={v[15298]} x5647={v[5647]} x34606={v[34606]} "
          f"| bad_checks={len(b)} failing_eqs={len(f)} score={L.NEQ-len(f)}")
    print(f"    checks: {b}")
