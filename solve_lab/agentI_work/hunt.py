#!/usr/bin/env python3
"""Prioritised hunt for the eq8680 compensator.

Candidates in order of how little they import.  Each is added to the
deliverable's support + a23434 and evaluated with exact branch-and-bound.
"""
import collections, itertools, json, os, sys, time
from eq8680 import build, minfail_bnb, M, wit, av, v2a, BASE, RES
from intsolve import solve_int

HERE = os.path.dirname(os.path.abspath(__file__))
E_BASE = sorted({e for a in BASE for e, _ in M.atom_eqs[a]})
SB = set(E_BASE)

CANDS = [
    ('X19964', (1631,)),          # net -1  -- the cheapest possible compensator
    ('X4432', (2427, 22331, 33706)),
    ('X6947', (23435,)),
    ('X33168', (23437,)),
    ('X10422', (11772,)),
    ('X11099', (11774,)),
    ('X22526', (11776,)),
    ('X34868', (11778,)),
    ('X950', (20290,)),
    ('X15120', (20292,)),
    ('X35531', (20294,)),
    ('X18253', (20292, 20293)),
    ('X37720', (20294, 20295)),
    ('X23822', (11776, 11777)),
    ('X7945', (11778, 11779)),
    ('X9629', (20290, 20291)),
    ('X37413', (11774, 11775)),
    ('X37254', (11775, 11917)),
    ('X35619', (23437, 23438, 35830)),
    ('X23642', (11772, 11773, 34814)),
    ('X30108', (11773, 34813, 36618)),
    ('X15324', (11775, 22331, 33707)),
]


def show(tag, add, tlimit):
    SUP = BASE + [a for a in add if a not in BASE]
    r = build(SUP)
    if r is None:
        print(f"  {tag} add={add}: non-linear knob -> skipped", flush=True)
        return None
    SUPs, knobs, E, base, Mat = r
    new = [e for e in E if e not in SB]
    t0 = time.time()
    mf, forced, nact, nodes = minfail_bnb(E, base, Mat, budget=6, tlimit=tlimit)
    if mf is None:
        txt = "minfail > 6  (cannot beat 39,026)"
    elif mf == 'timeout':
        txt = f"TIMEOUT ({nodes} nodes)"
    else:
        txt = f"minfail = {mf}"
    star = '   *** BEATS 39,026 ***' if isinstance(mf, int) and mf < 7 else ''
    print(f"  {tag:8s} add={str(add):28s} |E|={len(E):3d} (+{len(new):2d} new) "
          f"knobs={len(knobs):2d} forced={forced} active={nact} "
          f"nodes={nodes:6d} {int(time.time()-t0):4d}s  {txt}{star}", flush=True)
    return {'tag': tag, 'add': list(add), 'nE': len(E), 'new': len(new),
            'knobs': len(knobs), 'minfail': mf, 'forced': forced}


def main():
    tl = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    print(f"base support {BASE}: |E|={len(E_BASE)}", flush=True)
    out = []
    for tag, add in CANDS:
        r = show(tag, add, tl)
        if r:
            out.append(r)
        json.dump(out, open(os.path.join(HERE, 'hunt_result.json'), 'w'), indent=1)
    # best pairs of the cheapest single-atom groups
    print("\npairs of the cheapest groups:", flush=True)
    cheap = [c for c in CANDS[:12]]
    for (t1, a1), (t2, a2) in itertools.combinations(cheap, 2):
        r = show(f"{t1}+{t2}", tuple(sorted(set(a1) | set(a2))), tl)
        if r:
            out.append(r)
        json.dump(out, open(os.path.join(HERE, 'hunt_result.json'), 'w'), indent=1)
    print("\ndone", flush=True)


if __name__ == '__main__':
    main()
