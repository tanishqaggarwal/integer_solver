#!/usr/bin/env python3
"""Sanity: for RANDOM bit settings, does forward-eval (Z) violate ONLY twist-family
atoms? If the violated set is always within {atoms touching x_9770/x_3183/x_18274/
x_17728}, then the search space is clean: twist match <=> full solve. Also count
how many atoms ever float."""
import json
from confluent_eval5 import build5, make_forward
from propagate import atom_vars

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]

def main():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]

    def viol(setbits):
        val = solve(list(bestval), list(setbits))
        vs = []
        for a, poly in enumerate(A):
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s: vs.append(a)
        return set(vs)

    everfloat = set()
    # deterministic pseudo-random subsets of various sizes/mixes
    tests = []
    for seed in range(8):
        S22 = [BITS22[i] for i in range(22) if ((i*7+seed*13) % 3 == 0)]
        S233 = [bits233[(i*97+seed*31) % len(bits233)] for i in range(seed*3 % 20 + 2)]
        tests.append(sorted(set(S22 + S233)))
    for S in tests:
        vs = viol(S)
        everfloat |= vs
        n22 = sum(1 for b in S if b in set(BITS22)); n233 = len(S)-n22
        print(f"|S|={len(S)} ({n22}A+{n233}B): {len(vs)} violated: {sorted(vs)[:12]}", flush=True)
    print(f"\nUnion of ever-violated atoms across tests: {sorted(everfloat)}")
    print(f"count = {len(everfloat)}")

if __name__ == '__main__':
    main()
