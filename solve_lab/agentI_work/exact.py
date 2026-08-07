#!/usr/bin/env python3
"""EXHAUSTIVE proof that minfail >= 7 for the tightest supports.

greedy only gives an UPPER bound on minfail.  Here we enumerate every set of
<= 6 sacrificed equations and test integer solvability of the rest, which settles
"can this support fail at most 6?" outright.
"""
import itertools, json, os, sys, time
from search import build, RES, M
from intsolve import solve_int

HERE = os.path.dirname(os.path.abspath(__file__))
CAND = [
    [],                                  # the deliverable's own support
    [23434],
    [23434, 23435, 23436],
    [23436, 23437],
    [23437, 23438],
    [23436, 23438],
    [23434, 36224, 36223],
    [23434, 23435],
    [23435, 23436],
]
out = []
for extra in CAND:
    r = build(RES + extra)
    if r is None:
        print("skip", extra); continue
    SUP, knobs, E, base, Mat = r
    n = len(E)
    t0 = time.time()
    hit = None
    for k in range(0, 7):
        for drop in itertools.combinations(range(n), k):
            keep = [i for i in range(n) if i not in drop]
            if solve_int([Mat[i] for i in keep], [-base[i] for i in keep]) is not None:
                hit = (k, [E[i] for i in drop]); break
        if hit:
            break
    res = {'extra': extra, 'nE': n, 'knobs': len(knobs),
           'minfail_le6': hit[0] if hit else None,
           'proved_ge7': hit is None,
           'secs': round(time.time() - t0, 1)}
    out.append(res)
    if hit:
        print(f"extra={extra} |E|={n} knobs={len(knobs)}  *** minfail={hit[0]} "
              f"sacrificing {hit[1]}  ({res['secs']}s)", flush=True)
    else:
        print(f"extra={extra} |E|={n} knobs={len(knobs)}  PROVED minfail >= 7 "
              f"(all C(n,<=6) subsets unsolvable)  ({res['secs']}s)", flush=True)
    json.dump(out, open(os.path.join(HERE, 'exact_result.json'), 'w'), indent=1)
print("\ndone")
