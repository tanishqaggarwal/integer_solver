#!/usr/bin/env python3
"""Test TRUE linearity (mod P) of x_9770,x_3183 in 22 bits and x_18274,x_17728
in 233 bits. Mod-P forward eval never leaves stale values (division invertible),
so it reveals the genuine algebraic structure. If linear mod several primes,
the wire is a linear form (subset-sum) over the bits -> lattice attack applies."""
import json, time
from confluent_eval5 import build5, make_forward
from propagate import NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
PRIMES = [(1 << 61) - 1, (1 << 31) - 1, 2305843009213693951 if False else (1<<44)-3]

def build_seq():
    A, kind, info, seq, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    return A, kind, info, seq, bestval

def main():
    t0 = time.time()
    A, kind, info, seq, bestval = build_seq()
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]

    for P in PRIMES:
        solve = make_forward(kind, info, seq, bestval, mod=P)
        bm = [x % P for x in bestval]
        base = solve(list(bm), [])
        b0 = {w: base[w] for w in (9770, 3183, 18274, 17728)}
        single = {}
        for b in control:
            val = solve(list(bm), [b])
            single[b] = {w: val[w] for w in (9770, 3183, 18274, 17728)}

        def lintest(watch, pool, sizes, seedmul):
            allok = True
            for si, sz in enumerate(sizes):
                idx = [((i*seedmul + si*7919 + 3) % len(pool)) for i in range(sz)]
                S = sorted(set(pool[i] for i in idx))
                val = solve(list(bm), S)
                actual = val[watch] % P
                pred = (b0[watch] + sum((single[b][watch] - b0[watch]) for b in S)) % P
                ok = (actual == pred)
                if not ok: allok = False
                print(f"    x_{watch} |S|={len(S)}: {'LIN' if ok else 'NONLIN'}", flush=True)
            return allok

        print(f"\n=== P={P} ===", flush=True)
        print("  x_9770 (22 bits):"); r1 = lintest(9770, BITS22, [3,7,12,19], 13)
        print("  x_3183 (22 bits):"); r2 = lintest(3183, BITS22, [3,7,12,19], 17)
        print("  x_18274 (233 bits):"); r3 = lintest(18274, bits233, [3,10,50,120,211], 101)
        print("  x_17728 (233 bits):"); r4 = lintest(17728, bits233, [3,10,50,120,211], 103)
        print(f"  SUMMARY P={P}: 9770={r1} 3183={r2} 18274={r3} 17728={r4}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
