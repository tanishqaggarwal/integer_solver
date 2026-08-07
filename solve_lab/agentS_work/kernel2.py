"""Complement to kernel.py.

kernel.py's limitation, stated plainly: motion along the affine kernel changes a20215 only by
multiples of p (that IS the p*Z^2 result), so the residual's mod-p class is INVARIANT along the
kernel.  The membership answer therefore cannot change unless the *measured knob set* changes at
the displaced configuration.  It is a test of structural stability, not of the congruence.

This script supplies the part that does move the class: at each BFS image configuration -- which
by construction sits at a DIFFERENT mod-p 5-tuple -- solve the other rows exactly to land on a
near-solution there, then test membership of the target residual in the reachable lattice.
That is lat3.analyse, which is the correct formulation (lat5.py demanded ALL rows including the
targets, which is strictly stronger and so a weaker probe of this question).

Knob set and base are re-measured per configuration; base configuration for each is
triple8_seed + that image point's selector settings.
"""
import sys, json, collections, pickle, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C, lat3
import harness as H, engine as E
P = C.P

seen = pickle.load(open('bfs_image.pkl', 'rb'))
items = sorted(seen.items(), key=lambda kv: len(kv[1]))
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 14
print("BFS image points available: %d ; testing %d" % (len(items), LIMIT), flush=True)

blocked = solved = infeas = 0
classes = collections.Counter()
for i, (k, a) in enumerate(items[:LIMIT]):
    seed = dict(C.BASE); seed.update(a)
    classes[(k[3], k[4])] += 1
    t0 = time.time()
    try:
        r = lat3.analyse(seed, 'img%d(|on|=%d)' % (i, len(a)))
    except Exception as e:
        print("  img%d analyse ERR %s" % (i, type(e).__name__), flush=True); continue
    if r is None:
        infeas += 1
    elif r[0] is not None:
        solved += 1
        print("  *** img%d MEMBERSHIP SUCCEEDED -> FULL SOLVE PATH ***" % i, flush=True)
        json.dump({str(x): str(int(z)) for x, z in seed.items()}, open('S_kernel2_hit_%d.json' % i, 'w'))
    else:
        blocked += 1
    print("  [img%d done %.0fs] running: blocked=%d solved=%d other-rows-infeasible=%d"
          % (i, time.time() - t0, blocked, solved, infeas), flush=True)

print("\n=== SUMMARY (kernel2) ===", flush=True)
print("distinct (a20215,a28647) mod-p classes covered: %d" % len(classes), flush=True)
print("blocked=%d  solved=%d  other-rows-infeasible=%d" % (blocked, solved, infeas), flush=True)
if solved:
    print("=> obstruction DISSOLVES at some image points -- it is cfg0-local.", flush=True)
else:
    print("=> obstruction blocked at every image point where the other rows were solvable.", flush=True)
