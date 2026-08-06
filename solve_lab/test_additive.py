#!/usr/bin/env python3
"""Decisive test: is x_9770(S) == sum_{b in S} x_9770({b}) (a clean linear form /
subset-sum over the bits, with the empty-set 'base' being a stale gate value)?
Test for x_9770,x_3183 over 22 bits and x_18274,x_17728 over 233 bits, on RANDOM
subsets (deterministic pseudo-random via index arithmetic, no RNG)."""
import json, time
from confluent_eval5 import build5, make_forward
from propagate import NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

def main():
    t0 = time.time()
    A, kind, info, seq, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]

    # single-bit values
    single = {}
    for b in control:
        val = solve(list(bestval), [b])
        single[b] = {9770: val[9770], 3183: val[3183], 18274: val[18274], 17728: val[17728]}
    print(f"got single-bit values ({time.time()-t0:.0f}s)", flush=True)

    def check(watch, pool, subsets):
        okcount = 0
        for S in subsets:
            val = solve(list(bestval), S)
            actual = val[watch]
            summ = sum(single[b][watch] for b in S)
            ok = (actual == summ)
            okcount += ok
            if not ok and okcount < 3:
                print(f"    x_{watch} |S|={len(S)}: MISMATCH actual-sum={actual-summ}", flush=True)
        print(f"  x_{watch}: additive-over-values holds on {okcount}/{len(subsets)} subsets", flush=True)
        return okcount == len(subsets)

    # deterministic pseudo-random subsets
    def subs(pool, sizes, seedmul):
        out = []
        for si, sz in enumerate(sizes):
            idx = [( (i*seedmul + si*7919 + 3) % len(pool)) for i in range(sz)]
            S = sorted(set(pool[i] for i in idx))
            out.append(S)
        return out

    print("22-bit side:")
    check(9770, BITS22, subs(BITS22, [2,3,5,8,11,15,19], 13))
    check(3183, BITS22, subs(BITS22, [2,3,5,8,11,15,19], 17))
    print("233-bit side:")
    check(18274, bits233, subs(bits233, [2,3,5,10,30,80,150,211], 101))
    check(17728, bits233, subs(bits233, [2,3,5,10,30,80,150,211], 103))
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
