"""Channel U=0, V=1 with c*d = 1:
     x_15298 = 0   -> FIRST CORE dead
     x_38170 = 0   -> group-2 mirror dead
     x_5647  = 1   -> x_37892 = x_24908 (control x_19750), x_13682 = x_14853 (FREE)
   Conditions left: a688, a1618, a29539, a26731, and group-1 mirror (x_3719, x_25118).
   Controls: x_19750 (a688), x_14853 (a1618), x_14515 (a29539), x_16742 (a26731),
             and {16441, 22917, 31339, 33708} SPARE for the mirror -> 4 for 2.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]


def pins_of(bit):
    return [(rec[0], rec[1]) for rec in LD.get(str(bit), [])]


def state(BITS, th=(), load=True):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in dict(th).items():
        v[k] = x
    fw.forward(v)
    if load:
        for b in BITS:
            for a, X in pins_of(b):
                if X in dict(th):
                    continue
                x = fw.solve_lin(a, X, v)
                if x is not None:
                    v[X] = x
        fw.forward(v)
    return v


if __name__ == '__main__':
    for BITS in [(438, 91), (438, 1203), (490, 91), (1530, 2441), (2081, 91)]:
        v = state(BITS)
        bad = fw.bad_checks(v)
        av = L.all_atom_values(v)
        f = L.failing_eqs(av)
        print(f"BITS={BITS}: U={v[7715]} V={v[34554]} x15298={v[15298]} x5647={v[5647]} "
              f"ab={v[38170]} cd={v[3896]}")
        print(f"   x3719%p==0:{v[3719]%P==0} x25118%p==0:{v[25118]%P==0} "
              f"| bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}")
        print(f"   bad={bad}", flush=True)
