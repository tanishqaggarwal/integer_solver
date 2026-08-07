"""Deficiency-directed search for INDEPENDENT test cases.

Criterion (mine, applied three times now and the one that matters): a configuration is an
independent test case iff its OTHER rows are solvable AND its POST-SOLVE target class is distinct.
Measuring the class before the re-solve is the trap; it is not used anywhere here.

Shape bookkeeping per configuration: knobs K, other-rows M, kernel dim d,
rank = K - d, deficiency = M - rank.  deficiency > 0 was infeasible 21/21 in the logs.

Pool: the BFS image points (img4 -- the one existence proof -- came from here).  Sharded by
WK/NW so several workers can run at once.  Emits one machine-readable RESULT line per config so
the shards can be combined without re-parsing prose.
"""
import os, sys, json, collections, pickle, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C, lat3
import harness as H, engine as E
P = C.P

WK = int(os.environ.get('WK', 0)); NW = int(os.environ.get('NW', 1))
START = int(os.environ.get('START', 0))
seen = pickle.load(open('bfs_image.pkl', 'rb'))
items = sorted(seen.items(), key=lambda kv: len(kv[1]))
print("pool: %d BFS image points; worker %d/%d from index %d" % (len(items), WK, NW, START), flush=True)

for i, (k, a) in enumerate(items):
    if i < START or i % NW != WK:
        continue
    seed = dict(C.BASE); seed.update(a)
    t0 = time.time()
    try:
        r = lat3.analyse(seed, 'img%d(|on|=%d)' % (i, len(a)))
    except Exception as e:
        print("RESULT %d ERR %s" % (i, type(e).__name__), flush=True); continue
    if r is None:
        # other rows infeasible (or lattice trivial) -> not an independent case
        print("RESULT %d other-rows-infeasible - - -" % i, flush=True)
        continue
    y = r[0]
    # recover the post-solve class directly rather than scraping
    knobs = r[3]; aff = r[4]; v0 = r[5]; bad0 = r[6]
    n0 = r[2]
    ns = dict(seed)
    for j, f in enumerate(knobs):
        if n0[j]:
            ns[f] = v0[f] + n0[j]
    try:
        v = E.forward(ns)
        nsx = {'v': v, '__builtins__': {}}
        cls = tuple(eval(H.acodes[t], nsx) % P for t in (20215, 28647))
        bad = E.badatoms(v); ff = E.eqfails(bad)
    except Exception:
        cls = None; ff = []
    status = 'SOLVED' if y is not None else 'blocked'
    print("RESULT %d %s %s %s score=%s (%.0fs)"
          % (i, status, cls[0] if cls else '-', cls[1] if cls else '-',
             39033 - len(ff) if cls else '-', time.time() - t0), flush=True)
    if y is not None:
        json.dump({str(x): str(int(z)) for x, z in ns.items()}, open('S_dir_hit_%d.json' % i, 'w'))
        print("*** SOLVED at img%d -- endgame condition DISSOLVES ***" % i, flush=True)
