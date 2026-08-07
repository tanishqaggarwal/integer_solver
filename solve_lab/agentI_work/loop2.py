#!/usr/bin/env python3
"""No-good loop built on the fast recorded run (fprun.run).

Blacklists non-boolean decision variables in each conflict cone and retries.
Persists the blacklist after every round.
"""
import os, sys, json, time, pickle, collections
import fprun
from model import Model
HERE = os.path.dirname(os.path.abspath(__file__))
BL = os.path.join(HERE, 'blacklist2.json')


def boolvars(M):
    bv = set()
    for a in range(M.na):
        q = M.polys[a]
        vs = set()
        for m in q:
            vs |= set(m)
        if len(vs) == 1 and max(len(m) for m in q) == 2:
            bv |= vs
    return bv


def main():
    pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    bl = set(json.load(open(BL))) if os.path.exists(BL) else set()
    Mtmp = Model()
    bv = boolvars(Mtmp)
    del Mtmp
    print("boolean vars:", len(bv), "blacklist:", len(bl), flush=True)
    hist = []
    for it in range(rounds):
        out, M, E = fprun.run(pol, seed, bl, verbose=False, tag=f'loop_{pol}')
        if out['status'] == 'ok':
            print(f"[{it}] *** COMPLETE MOD-P SOLUTION *** decisions={len(out['decisions'])}")
            pickle.dump(out, open(os.path.join(HERE, f'fp_complete_{pol}_{seed}.pkl'), 'wb'))
            return
        decs = set(out['decs'])
        new = sorted(v for v in decs if v not in bv and v not in bl)
        print(f"[{it}] conflict a{out['bad']} decs={len(decs)} new={len(new)} "
              f"known={sum(1 for x in out['val'] if x is not None)}", flush=True)
        hist.append((out['bad'], sorted(decs)))
        if not new:
            print("   STUCK.  cone decisions:", sorted(decs)[:30])
            print("   boolean decisions in cone:", sorted(v for v in decs if v in bv)[:30])
            json.dump({'stuck_at': out['bad'], 'decs': sorted(decs)},
                      open(os.path.join(HERE, 'stuck.json'), 'w'))
            break
        bl |= set(new)
        json.dump(sorted(bl), open(BL, 'w'))
    json.dump(sorted(bl), open(BL, 'w'))


if __name__ == '__main__':
    main()
