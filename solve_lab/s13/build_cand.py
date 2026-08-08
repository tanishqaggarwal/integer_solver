#!/usr/bin/env python3
"""
Build a full assignment from the clean-frame knob solution and verify it with
the EXACT integer checker.  This is the unconditional test (T5): whatever the
linear model claims, only checker.py decides.

Usage: python3 build_cand.py
"""
import os, sys, json, time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import fwd_frame as F
from reduce_tw import downstream

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977


def main():
    t0 = time.time()
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    sol = json.load(open(os.path.join(HERE, 'knob_solution.json')))
    print(f"knob solution: {len(sol)} knobs (mod p)")

    unknown, touched = downstream([8731, 9118])
    newv = list(v)
    changed = []
    for k, t in sol.items():
        i = int(k[2:])
        if newv[i] % P != t % P:
            changed.append((i, newv[i] % P, t % P))
        # keep the same k*p offset, change only the residue
        newv[i] = newv[i] - (newv[i] % P) + (t % P)
    print(f"knobs whose residue actually changes: {len(changed)}")
    for i, a, b in changed[:10]:
        print(f"   x{i}: {str(a)[:22]}... -> {str(b)[:22]}...")

    # clean-frame forward evaluation over the cone
    val = F.evaluate(newv, {}, unknown, touched, frozen={})
    for x, t in val.items():
        newv[x] = t

    out = os.path.join(HERE, 'cand_clean.json')
    json.dump({f'x_{i}': int(newv[i]) for i in range(len(newv)) if newv[i]},
              open(out, 'w'))
    print(f"\nwritten -> {out}")
    print(f"{time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
