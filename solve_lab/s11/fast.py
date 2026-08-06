import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
BITS = (542, 47, 438, 91)
NAMES = ['x3719', 'x25118', 'x25614', 'x34220', 'n-gap', 'm-gap']


def light(th):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in th.items():
        v[k] = x
    fw.forward(v)
    return v


def targets(v):
    return [v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P,
            (v[12186] - v[1308]) % P, (v[24908] - v[19083]) % P]


POOL = [14515, 19750, 5096, 21589, 33708, 31339, 28486, 29261, 26489, 2467, 19275, 28548,
        6250, 5460, 8363, 30060, 32184, 3271, 18944, 1962, 8971, 3473, 8060, 5616, 245, 19450,
        27156, 12871, 7994, 18288, 14806, 37589, 23342, 26510, 3812, 33900, 17983, 5799, 8431,
        36711, 11258, 32322]

if __name__ == '__main__':
    th = {}
    r = targets(light(th))
    t0 = time.time()
    rows = {i: [] for i in range(6)}
    for c in POOL:
        t2 = dict(th)
        t2[c] = 1
        r1 = targets(light(t2))
        for i in range(6):
            if (r1[i] - r[i]) % P:
                rows[i].append(c)
    print(f"({time.time()-t0:.0f}s)")
    for i in range(6):
        print(f"  {NAMES[i]:8s}: {rows[i]}")
