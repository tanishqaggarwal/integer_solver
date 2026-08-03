#!/usr/bin/env python3
"""Test the exact multilinear degree (mod P) of the numerator-core wires x_15690,
x_26870 (x_6773 = x_15690 - x_26870) and x_14494, x_26234 (x_17233 side), plus the
shared denominator x_8821, over the 233 bits. If x_15690/x_14494 are degree<=2, the
twist equation (for a FIXED 22-side target N) is quadratic in the 233 bits ->
linearization/lattice attack becomes possible."""
import json, time, itertools
from confluent_eval5 import build5, make_forward
P = (1 << 61) - 1
BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])
WIRES = [15690, 26870, 14494, 26234, 8821, 6773, 17233]

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]
    probe = bits233[:12]

    cache = {}
    def ev(S):
        key = tuple(sorted(S))
        if key in cache: return cache[key]
        val = solveP(list(bm), list(S))
        r = {w: val[w] for w in WIRES}
        cache[key] = r; return r

    f0 = ev([])
    fs = {b: ev([b]) for b in probe}
    pair = {}
    c2 = {w: 0 for w in WIRES}
    for b, c in itertools.combinations(probe, 2):
        v = ev([b, c]); pair[(b, c)] = v
        for w in WIRES:
            if (v[w] - fs[b][w] - fs[c][w] + f0[w]) % P: c2[w] += 1
    npair = len(list(itertools.combinations(probe, 2)))
    print(f"deg2 nonzero /{npair}: " + ", ".join(f"x_{w}={c2[w]}" for w in WIRES), flush=True)
    c3 = {w: 0 for w in WIRES}
    ntri = 0
    for b, c, d in itertools.combinations(probe, 3):
        ntri += 1
        vv = ev([b, c, d])
        for w in WIRES:
            m3 = (vv[w] - pair[(b,c)][w] - pair[(b,d)][w] - pair[(c,d)][w]
                  + fs[b][w] + fs[c][w] + fs[d][w] - f0[w]) % P
            if m3: c3[w] += 1
    print(f"deg3 nonzero /{ntri}: " + ", ".join(f"x_{w}={c3[w]}" for w in WIRES), flush=True)
    c4 = {w: 0 for w in WIRES}; nq = 0
    for Q in itertools.combinations(probe[:9], 4):
        nq += 1
        s = {w: 0 for w in WIRES}
        for r in range(5):
            for sub in itertools.combinations(Q, r):
                vv = ev(list(sub))
                for w in WIRES: s[w] = (s[w] + (-1)**(4-r)*vv[w]) % P
        for w in WIRES:
            if s[w]: c4[w] += 1
    print(f"deg4 nonzero /{nq}: " + ", ".join(f"x_{w}={c4[w]}" for w in WIRES), flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
