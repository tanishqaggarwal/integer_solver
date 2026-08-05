#!/usr/bin/env python3
"""The 22-side produces only 174 distinct x_9770 values. Get the EXACT Z values and
analyze structure: arithmetic progression? small set of residue sums? This reveals
the codeword structure the 233-side must also hit."""
import json, time
import numpy as np
from math import gcd
from functools import reduce
from confluent_eval5 import build5, make_forward
from propagate import NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
p = 2147483647

def main():
    t0 = time.time()
    a9 = np.load('tab22_9770_2147483647.npy'); a3 = np.load('tab22_3183_2147483647.npy')
    # find one representative A (index) per distinct x_9770 mod p
    rep = {}
    for A in range(len(a9)):
        v = int(a9[A])
        if v not in rep: rep[v] = A
    print(f"distinct x_9770 mod p: {len(rep)} ({time.time()-t0:.0f}s)", flush=True)

    Aatoms, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)

    # exact Z x_9770 for each representative
    vals9770 = []
    vals3183 = []
    for vmod, Aidx in rep.items():
        Abits = [BITS22[i] for i in range(22) if (Aidx >> i) & 1]
        val = solve(list(bestval), Abits)
        vals9770.append(val[9770]); vals3183.append(val[3183])
    vals9770 = sorted(set(vals9770)); vals3183 = sorted(set(vals3183))
    print(f"exact distinct x_9770: {len(vals9770)}; x_3183: {len(vals3183)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"x_9770 range: min={vals9770[0]}\n              max={vals9770[-1]}", flush=True)
    # is it an arithmetic progression? check consecutive differences
    diffs = sorted(set(vals9770[i+1]-vals9770[i] for i in range(len(vals9770)-1)))
    print(f"x_9770 sorted-consecutive diffs: {len(diffs)} distinct; smallest few: {diffs[:5]}", flush=True)
    g = reduce(gcd, [abs(v) for v in vals9770 if v] + [vals9770[i+1]-vals9770[i] for i in range(len(vals9770)-1)])
    print(f"gcd of all x_9770 values & diffs: {g}", flush=True)
    if g > 1:
        base = vals9770[0]
        ks = sorted(set((v-base)//g for v in vals9770))
        print(f"  as base + k*g: k values span {ks[0]}..{ks[-1]}, {len(ks)} distinct; contiguous={ks==list(range(ks[0],ks[-1]+1))}", flush=True)
        print(f"  base={base}, step g={g}", flush=True)
    # save exact codewords
    json.dump({'x9770': [str(v) for v in vals9770], 'x3183': [str(v) for v in vals3183]}, open('codewords.json','w'))
    print(f"wrote codewords.json ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
