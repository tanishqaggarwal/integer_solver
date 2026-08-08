#!/usr/bin/env python3
"""deep_frontier.py -- pin the TINY success rates at the wall edge with many
independent tabu restarts, so F (any nonzero rate) and F (>1%) are settled and
not an artifact of a lucky 1-in-2.  Also records best-E, and the same for the
one-operand-known (clamp a) sub-instance."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import fbuild as FB          # sets sys.path to reach solvers
import solvers as S

TIME_CAP = 75.0            # seconds per (s, mode) cell
ITERS = 150000            # tabu high effort
OUT = os.path.join(HERE, 'deep_frontier.json')


def deep(s, clamp_which=None):
    mm = FB.baseline_modmul(s, mode='wallace')
    x0, clamp = (None, ())
    if clamp_which:
        x0, clamp = FB.clamp_operand(mm, which=clamp_which,
                                     rng=np.random.default_rng(7))
    hits = 0; best = 1e18; n = 0; t0 = time.time()
    for r in range(200000):
        e, _ = S.tabu(mm['ising'], iters=ITERS, seed=100000 + r, x0=x0, clamp=clamp)
        best = min(best, e); hits += (e == 0); n += 1
        if time.time() - t0 > TIME_CAP:
            break
    return dict(s=s, n_vars=mm['Q'].n, hits=int(hits), n=int(n),
                rate=hits / n, best=float(best), elapsed=round(time.time()-t0, 1),
                clamp=clamp_which or 'free')


if __name__ == '__main__':
    res = {}
    for tag, cw in [('free', None), ('clampA', 'a')]:
        for s in [6, 7, 8, 9, 10]:
            d = deep(s, clamp_which=cw)
            res[f"{tag}|{s}"] = d
            print(f"[{tag}] s={s:2d} n={d['n_vars']:5d} tabu {d['hits']:3d}/{d['n']:<5d} "
                  f"rate={d['rate']*100:6.3f}% best=E{d['best']:.0f} ({d['elapsed']}s)", flush=True)
            json.dump(res, open(OUT, 'w'), indent=1)
    print("DEEP DONE")
