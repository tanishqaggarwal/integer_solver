"""Isolate the load-pin cascade: set each bit's data inputs, see exactly what breaks."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
LD = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'loads.json')))['loads']


def pins_of(bit):
    return [(rec[0], rec[1]) for rec in LD[str(bit)]]


for BITS in [(542, 438), (24601, 438), (47, 438)]:
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    fw.forward(v)
    before = set(fw.bad_checks(v))
    print(f"=== BITS={BITS}  bad before loading = {len(before)} {sorted(before)}", flush=True)
    for b in BITS:
        for a, X in pins_of(b):
            x = fw.solve_lin(a, X, v)
            ok = x is not None
            if ok:
                v[X] = x
            fw.forward(v)
            print(f"    bit {b}: pin a{a} -> set x{X} {'OK' if ok else 'FAIL'} ; "
                  f"bad now {len(fw.bad_checks(v))}", flush=True)
    after = fw.bad_checks(v)
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    print(f"    after loading: bad={len(after)} failing={len(f)} score={L.NEQ-len(f)}")
    print(f"      new bad: {[a for a in after if a not in before]}")
    print(f"      resolved: {[a for a in before if a not in set(after)]}", flush=True)
