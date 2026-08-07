"""Inspect the 39,026 checkpoint in the session-11 structural coordinates."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
LAB = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
d = json.load(open(os.path.join(LAB, 'best', 'new_instance_partial_39026.json')))
v = [0] * L.NVARS
for k, x in d.items():
    v[int(k[2:]) if k.startswith('x_') else int(k)] = int(x)
av = L.all_atom_values(v)
nz = [a for a, x in enumerate(av) if x != 0]
f = L.failing_eqs(av)
print(f"checkpoint: nonzero atoms={len(nz)} {nz}   failing_eqs={len(f)} score={L.NEQ-len(f)}")
bad = [a for a in nz if L.atom_out.get(a) is None]
gat = [a for a in nz if L.atom_out.get(a) is not None]
print(f"  of those: CHECK atoms={bad}  GATE atoms (broken gates)={gat}")
print()
print("STRUCTURAL COORDINATES AT THE CHECKPOINT")
for name, u in [('x_8599 (a)', 8599), ('x_21839 (b)', 21839), ('x_7304 (c)', 7304), ('x_25956 (d)', 25956),
                ('U = x_7715', 7715), ('V = x_34554', 34554),
                ('x_15298 = UV', 15298), ('x_5647 = (1-U)V', 5647), ('x_34606 = U(1-V)', 34606),
                ('x_38170 = a*b', 38170), ('x_3896 = c*d', 3896),
                ('x_4287', 4287), ('x_2081', 2081), ('x_13195', 13195)]:
    print(f"  {name:20s} = {v[u]}")
print()
print("  x_37892 =", str(v[37892])[:50], " x_13682 =", str(v[13682])[:50])
print("  x_12186 - x_1308  =", (v[12186]-v[1308]) % P == 0, " (mod p zero?)")
print("  x_24908 - x_19083 =", (v[24908]-v[19083]) % P == 0)
for q in (3719, 25118, 25614, 34220):
    print(f"  x_{q} mod p == 0 : {v[q] % P == 0}")
print()
print("how many boolean message bits are ON at the checkpoint:")
ld = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'loads.json')))['loads']
on = [int(b) for b in ld if v[int(b)] != 0]
print(f"  {len(on)} of {len(ld)} : {on[:30]}")
