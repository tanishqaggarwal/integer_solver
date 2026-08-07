"""Channel U=1, V=0 with a*b = 1:
     x_15298 = 0  -> FIRST CORE structurally dead
     x_3896  = 0  -> group-1 mirror dead
     x_34606 = 1  -> x_37892 = x_16742 (free), x_13682 = x_12186 (controllable via x_5096)
   Remaining: group-2 mirror (2 conditions) with 4 spare controls, + a26731, + load pins.
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


def state(BITS, th):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in th.items():
        v[k] = x
    fw.forward(v)
    # load the active bits' constants
    for b in BITS:
        for a, X in pins_of(b):
            if X in th:
                continue
            x = fw.solve_lin(a, X, v)
            if x is not None:
                v[X] = x
    fw.forward(v)
    return v


def report(BITS, th, tag=''):
    v = state(BITS, th)
    bad = fw.bad_checks(v)
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    print(f"  {tag} U={v[7715]} V={v[34554]} x15298={v[15298]} x34606={v[34606]} ab={v[38170]} cd={v[3896]}")
    print(f"     x25614%p==0:{v[25614]%P==0} x34220%p==0:{v[34220]%P==0} "
          f"a26731gap={(v[16742]-v[19083])%P==0} a29539gap={(v[14853]-v[1308])%P==0}")
    print(f"     bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}  {bad}", flush=True)
    return v, bad, f


if __name__ == '__main__':
    for BITS in [(542, 47), (542, 1502), (1685, 47), (853, 112)]:
        print(f"=== BITS={BITS}")
        report(BITS, {}, 'raw   ')
