#!/usr/bin/env python3
"""agent V -- V6.  THE GROWTH CURVE: how many handle-less atoms end up nonzero, versus |S|?

The question this decides: is the handle-less population a nuisance that stays small at large |S|,
or does it switch on at a threshold and become a wall?

Measured at the CONSTRUCTOR's greedy fixpoint -- the state closeS4's closure loop starts from --
because that is where the population is set: the closure loop's global guard only ever accepts a
shift that strictly DECREASES the total nonzero count, so it cannot create handle-less nonzeros,
only remove them.  The fixpoint count is therefore an upper bound on what any run of that loop can
leave behind, and it is ~5 s per configuration instead of ~200 s.

Scope note the lab keeps needing: FOUR ON-sets per size, not one.  Seed 7 is L's own convention so
the row is comparable with T's numbers; seeds 101/202/303 measure whether the size or the draw is
doing the work.

Usage: python3 v_hcurve.py [outfile-tag]
"""
import os, sys, json, time, random, collections
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p = B.E, B.SL, B.p
atomh = B.GL['atomh']; M = B.GL['M']
HLS = set(a for a in E.res if len(atomh[a]) == 0)

SIZES = [1, 2, 4, 8, 12, 16, 17, 20, 24, 32, 40, 48, 64, 96, 128, 176, 224, 256]
SEEDS = [7, 101, 202, 303]

if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'curve'
    open(os.path.join(V, 'run_%s.pid' % tag), 'w').write('%d\n' % os.getpid())
    print('pid %d ; %d handle-less atoms of %d' % (os.getpid(), len(HLS), len(E.res)), flush=True)
    print('%-6s %-8s %-10s %-10s %-10s  %s' % ('|S|', 'seed', 'globalNZ', 'handleless', 'c>1 viol',
                                               'which handle-less'), flush=True)
    rows = []
    seen = collections.Counter()
    t0 = time.time()
    for n in SIZES:
        for sd in SEEDS:
            if n == 2 and sd == 7:
                S = B.onset(2)
            elif n > len(M['live']):
                continue
            else:
                S = random.Random(sd).sample(M['live'], n) if sd != 7 else B.onset(n)
            vv = B.greedy_init(S)
            r = E.run(vv)
            nz = [E.res[i] for i, x in enumerate(r) if x]
            hl = [a for a in nz if a in HLS]
            viol = len(B.violated(vv, r))
            for a in hl:
                seen[a] += 1
            rows.append({'n': n, 'seed': sd, 'global_nz': len(nz), 'handleless_nz': len(hl),
                         'viol': viol, 'which': hl})
            print('%-6d %-8d %-10d %-10d %-10d  %s'
                  % (n, sd, len(nz), len(hl), viol, [x[:40] for x in hl]), flush=True)
            json.dump({'rows': rows, 'seen': dict(seen)},
                      open(os.path.join(V, 'v_hcurve_%s.json' % tag), 'w'), indent=1)
    print('\ndistinct handle-less atoms ever nonzero: %d' % len(seen), flush=True)
    for a, k in seen.most_common():
        print('   %3d occurrence(s)  %s' % (k, a[:100]), flush=True)
    print('\nwall %.1f s' % (time.time()-t0), flush=True)
