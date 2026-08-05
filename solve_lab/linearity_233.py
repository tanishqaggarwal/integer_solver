#!/usr/bin/env python3
"""Is x_18274 (and x_17728) LINEAR in the 233 bits under forward-eval? If
x_18274 = base + sum_i b_i*delta_i exactly, then 'x_18274 = target' is an integer
subset-sum sum_i b_i*delta_i = target-base -> attack with lattice/LLL. Extract the
per-bit increments delta_i, check linearity on random combos, and report the gcd /
bit-size structure of the deltas (which controls lattice density)."""
import json, time
from math import gcd
from functools import reduce
from confluent_eval5 import build5, make_forward

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    b233 = [b for b in control if b not in BITS22]
    b22 = [b for b in control if b in BITS22]

    base = solve(list(bestval), [])
    for W in (18274, 17728, 9770, 3183):
        side = b233 if W in (18274, 17728) else b22
        delta = {}
        for b in side:
            v = solve(list(bestval), [b])
            delta[b] = v[W] - base[W]
        # linearity test: random subsets, predicted vs actual
        st = 999
        def rnd():
            nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
        nlin = 0; ntest = 40
        for _ in range(ntest):
            k = 1 + rnd() % min(10, len(side))
            S = sorted(set(side[rnd() % len(side)] for _ in range(k)))
            v = solve(list(bestval), S)
            pred = base[W] + sum(delta[b] for b in S)
            if v[W] == pred: nlin += 1
        nz = [d for d in delta.values() if d != 0]
        g = reduce(gcd, [abs(d) for d in nz]) if nz else 0
        bits = [abs(d).bit_length() for d in nz]
        print(f"x_{W}: linear on {nlin}/{ntest} random subsets; {len(nz)}/{len(side)} bits move it", flush=True)
        print(f"   base={base[W]} ({base[W].bit_length()} bits)", flush=True)
        print(f"   gcd(deltas)={g} ({g.bit_length()} bits); delta bit-sizes min/max={min(bits) if bits else 0}/{max(bits) if bits else 0}", flush=True)
        if g and nz:
            red = [d//g for d in nz]
            rb = [abs(r).bit_length() for r in red]
            print(f"   deltas/gcd bit-sizes min/max = {min(rb)}/{max(rb)}; count={len(red)}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
