#!/usr/bin/env python3
"""With the slack gate ON (x_12779=1 via a 22-side bit), test whether the bridge
wire x_24026 and the twist gap are LINEAR in the 233 bits (mod P). If linear, the
bridge equation x_24026 = x_18274 - x_35186 is a subset-sum -> lattice-solvable.
Test several x_12779=1 base bits."""
import json, time
from confluent_eval5 import build5, make_forward
P = (1 << 61) - 1
BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

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

    WATCH = [24026, 12520, 29798, 9770, 18274, 3368, 12779]
    for basebit in [19520, 2795, 1858]:
        base = solveP(list(bm), [basebit])
        print(f"\n== base bit {basebit}: x_12779={base[12779]}, x_24026={base[24026]}, x_9770={base[9770]} ==", flush=True)
        single = {}
        for b in bits233:
            single[b] = solveP(list(bm), [basebit, b])
        # linearity test of x_24026 and x_9770 over 233 bits (added to base)
        def lintest(w, sizes, seedmul):
            out = []
            for si, sz in enumerate(sizes):
                S = sorted(set(bits233[(i*seedmul+si*131+3) % len(bits233)] for i in range(sz)))
                val = solveP(list(bm), [basebit]+S)
                actual = val[w]
                pred = (base[w] + sum((single[b][w]-base[w]) for b in S)) % P
                out.append('LIN' if actual == pred else 'NL')
            return out
        for w in (24026, 9770, 18274):
            r = lintest(w, [3, 10, 40, 120, 233], 91)
            print(f"   x_{w}: {r}", flush=True)
        # how many 233-bits activate x_24026 (make it nonzero) from this base?
        act = [b for b in bits233 if single[b][24026] != base[24026]]
        print(f"   233-bits moving x_24026: {len(act)}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
