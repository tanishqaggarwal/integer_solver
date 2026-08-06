"""The cheap channel: U=V=1 with BOTH mirror gates off (a*b = 0, c*d = 0)."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
# a-bits set x_8599, b-bits x_21839, c-bits x_7304, d-bits x_25956
CFG = {
    'a=1,c=1 (542,438)': (542, 438),
    'a=1,d=1 (542,91)': (542, 91),
    'b=1,c=1 (47,438)': (47, 438),
    'b=1,d=1 (47,91)': (47, 91),
    'checkpoint-like b,c (24601?,438)': (24601, 438),
}
for name, BITS in CFG.items():
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    fw.forward(v)
    gaps = ((v[12186] - v[1308]) % P, (v[24908] - v[19083]) % P)
    core = (v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P)
    bad = fw.bad_checks(v)
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    print(f"{name}: U={v[7715]} V={v[34554]} x15298={v[15298]} ab={v[38170]} cd={v[3896]}")
    print(f"    gaps(zero?)={[g==0 for g in gaps]}  cores(zero?)={[c==0 for c in core]}")
    print(f"    bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}  {bad}")
