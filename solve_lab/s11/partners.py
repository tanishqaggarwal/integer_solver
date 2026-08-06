import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, leaves
P = L.P
BITS = (542, 47, 438, 91)
# linking checks:  a21050: x16441 ~ x4920 ; a34580: x33708 ~ x10170 ; a33796: x31339 ~ x6858
PAIRS = [(16441, 4920), (33708, 10170), (31339, 6858)]
v = [0] * L.NVARS
for b in BITS:
    v[b] = 1
fw.forward(v)
for ctl, partner in PAIRS:
    fr, _ = leaves.cone_free(partner)
    live = []
    for u in sorted(fr):
        old = v[u]
        v[u] = old + 1
        fw.forward(v)
        d = v[partner]
        v[u] = old
        fw.forward(v)
        if (d - v[partner]) % P:
            live.append(u)
    print(f"x{partner} (partner of x{ctl}): cone free={len(fr)}, MOD-P LIVE controls = {live}", flush=True)
