#!/usr/bin/env python3
"""Fit a generative model that exactly reproduces a reference assignment:
iteratively demote mis-derived variables to 'free' until forward evaluation
reproduces the reference everywhere."""
import os, pickle, json, sys, time
import jengine as E

HERE = os.path.dirname(os.path.abspath(__file__))
ref = E.load(sys.argv[1] if len(sys.argv) > 1 else
             os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
definer = E.build_definer()
order, cyc = E.topo(definer)
assert not cyc
ev = E.make_eval(definer)
free = set(range(E.NV)) - set(definer)
t0 = time.time()
for it in range(60):
    val = list(ref)
    # zero out derived vars to be safe
    bad = E.forward(val, order, ev, definer, free)
    diff = [i for i in range(E.NV) if val[i] != ref[i]]
    print(f"iter {it}: free={len(free)} inexact={len(bad)} mismatch={len(diff)} ({time.time()-t0:.0f}s)")
    if not diff and not bad:
        break
    free |= set(diff) | set(bad)
s, fails, av = E.score(val)
print("score of reproduced state:", s, fails)
print("final free count:", len(free))
pickle.dump({'definer': definer, 'order': order, 'free': free},
            open(os.path.join(HERE, 'jfit.pkl'), 'wb'))
