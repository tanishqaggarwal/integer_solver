#!/usr/bin/env python3
"""For each control bit, toggle it alone (from all-0) and count how many atoms the
Z forward-eval floats. 'Internally safe' = only twist-family atoms float (i.e. the
bit keeps every 233/22-side division wire exact). This tells us how constrained the
integer-consistent set is."""
import json
from confluent_eval5 import build5, make_forward

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])
TWIST = {1817, 30378, 40782, 44271}

def main():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))

    def nfloat(setbits):
        val = solve(list(bestval), list(setbits))
        extra = 0; tw = 0
        for a, poly in enumerate(A):
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s:
                if a in TWIST: tw += 1
                else: extra += 1
        return tw, extra

    safe22 = []; safe233 = []; unsafe = []
    for b in control:
        tw, extra = nfloat([b])
        side = '22' if b in BITS22 else '233'
        if extra == 0:
            (safe22 if b in BITS22 else safe233).append(b)
        else:
            unsafe.append((b, side, extra))
    print(f"internally-SAFE single bits: 22-side={len(safe22)}/22, 233-side={len(safe233)}/233")
    print(f"  safe22: {sorted(safe22)}")
    print(f"  safe233 ({len(safe233)}): {sorted(safe233)[:60]}")
    print(f"UNSAFE single bits (break internal divisions): {len(unsafe)}")
    # distribution of how many extra atoms unsafe bits break
    from collections import Counter
    dist = Counter(e for _, _, e in unsafe)
    print(f"  extra-atom-count distribution among unsafe: {dict(sorted(dist.items()))}")
    json.dump({'safe22': sorted(safe22), 'safe233': sorted(safe233),
               'unsafe': [[b, s, e] for b, s, e in unsafe]}, open('safe_bits.json', 'w'))
    print("wrote safe_bits.json")

if __name__ == '__main__':
    main()
