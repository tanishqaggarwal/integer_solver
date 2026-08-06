#!/usr/bin/env python3
"""MOD-P (artifact-free) how many atoms float for random B? If only the twist (4),
the 233 bits are UNconstrained mod-P (claw-find over full 2^233 = trapdoor). If more
float, the combination checks constrain B => consistent set may be small/enumerable."""
import json
from confluent_eval5 import build5, make_forward
P = (1 << 61) - 1
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
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]

    def floatP(setbits):
        val = solveP(list(bm), list(setbits))
        tw = ex = 0; exatoms = []
        for a, poly in enumerate(A):
            s = 0
            for m, c in poly.items():
                t = c % P
                for x in m: t = (t*val[x]) % P
                s = (s+t) % P
            if s % P:
                if a in TWIST: tw += 1
                else:
                    ex += 1
                    if len(exatoms) < 15: exatoms.append(a)
        return tw, ex, exatoms

    print(f"all-0: {floatP([])[:2]}", flush=True)
    for seed in range(6):
        S = sorted(set(bits233[(i*97+seed*31) % len(bits233)] for i in range(seed*4+3)))
        tw, ex, exa = floatP(S)
        print(f"|B|={len(S)}: twist={tw}, EXTRA={ex} {exa}", flush=True)
    # single bits
    xs = 0
    for b in bits233[:20]:
        tw, ex, _ = floatP([b])
        if ex: xs += 1
    print(f"single 233-bits (first 20) with extra mod-P floats: {xs}/20", flush=True)

if __name__ == '__main__':
    main()
