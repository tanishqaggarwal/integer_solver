#!/usr/bin/env python3
"""Does the CHOICE of dependent variable in the cascade closure matter?

The closure consumes one variable as a dependent for every atom it absorbs.  The
greedy version picks the first usable candidate.  Here the choice is randomised
many times to test whether any closure yields MORE than 7 knobs on the
deliverable's support -- i.e. whether the choice, not the rule, was the limit.
"""
import collections, random, sys, json, os
from cascade import M, wit, av, v2a, RES, unit_solvable

HERE = os.path.dirname(os.path.abspath(__file__))


def close_rand(SUP, rng, cap=4000):
    INSIDE = set(SUP)
    R = set(); dep = {}
    cnt_out = {y: sum(1 for a in v2a[y] if a not in INSIDE) for y in v2a}
    work = list(range(M.na))
    rng.shuffle(work)
    work = collections.deque(work)
    inq = bytearray(b'\x01') * M.na
    while work and len(R) < cap:
        b = work.popleft(); inq[b] = 0
        if b in INSIDE:
            continue
        cands = [y for y in M.avars[b]
                 if y not in dep and cnt_out[y] == 1 and b in v2a[y]
                 and unit_solvable(b, y)]
        if not cands:
            continue
        pick = rng.choice(cands)
        R.add(b); INSIDE.add(b); dep[pick] = b
        for y in M.avars[b]:
            cnt_out[y] -= 1
            for a2 in v2a[y]:
                if a2 not in INSIDE and not inq[a2]:
                    inq[a2] = 1; work.append(a2)
    knobs = sorted(y for y in v2a if y not in dep and cnt_out[y] == 0)
    return R, dep, knobs


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    best = None
    seen = collections.Counter()
    for s in range(N):
        rng = random.Random(s)
        R, dep, knobs = close_rand(RES, rng)
        seen[(len(R), len(knobs), tuple(knobs))] += 1
        if best is None or len(knobs) > best[0]:
            best = (len(knobs), len(R), knobs)
        print(f"  seed {s:3d}: absorbed {len(R):5d} atoms, {len(dep):5d} dependents, "
              f"KNOBS = {len(knobs)}  {knobs}", flush=True)
    print(f"\nmax knobs over {N} randomised closures: {best[0]}")
    print(f"distinct (|R|, #knobs, knobset) outcomes: {len(seen)}")
    json.dump({'max_knobs': best[0], 'knobs': best[2], 'trials': N},
              open(os.path.join(HERE, 'cascade_rand.json'), 'w'))


if __name__ == '__main__':
    main()
