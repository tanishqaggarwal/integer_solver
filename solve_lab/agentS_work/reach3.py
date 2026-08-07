"""sec 3 is refuted (reach2.py).  sec 2's JOINT condition is the statement that survives:
   x_18956 moves a20215 by 1 but pays 8863713 into a747, and a747's only handle steps by p, so
   keeping a747 satisfied forces 8863713*n = 0 (mod p), i.e. n = 0 (mod p) -- a20215 moves by
   multiples of p ONLY once the other rows must stay satisfied.  That is the pZ^2 image.

   Open question this settles: is that joint obstruction configuration-INDEPENDENT, or is it also
   just a cfg0-local fact?  Test it at random configurations OUTSIDE the BFS closure.
"""
import sys, json, collections, pickle, random, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C, lat2
import harness as H, engine as E, fast, sparse
P = C.P
BOOLS = sorted(f for f in C.cluster_cone() if C.isbool(f))
random.seed(4242)
cfgs = []
for sz in [17, 64, 128, 200]:
    S = set(random.sample(BOOLS, sz))
    a = {f: (1 if f in S else 0) for f in BOOLS}
    a[30163] = random.choice([0, 1]); a[11559] = random.choice([0, 1])
    cfgs.append(('rand|S|=%d' % sz, a))

for tag, extra in cfgs:
    seed = dict(C.BASE); seed.update(extra)
    t0 = time.time()
    try:
        v0, bad0, aff, atoms, hs = lat2.system(seed)
    except Exception as e:
        print("[%s] system ERR %s" % (tag, type(e).__name__), flush=True); continue
    nh = [x for x in atoms if x not in hs]
    knobs = sorted(aff)
    rows = [{f: aff[f][x] for f in knobs if x in aff[f]} for x in atoms]
    rhs = [-bad0.get(x, 0) for x in atoms]
    sol, msg, _ = sparse.solve_sparse(rows, rhs, names=atoms, verbose=False,
                                      maxcore=600, maxcorebits=400_000)
    print("[%s] bad=%d atoms=%d knobs=%d nohandle=%s -> %s  (%.0fs)"
          % (tag, len(bad0), len(atoms), len(knobs), nh,
             'FEASIBLE!!!' if sol is not None else msg, time.time() - t0), flush=True)
    if sol is not None:
        ns = dict(seed)
        for f, dv in sol.items():
            if dv: ns[f] = v0[f] + dv
        v = E.forward(ns); bad = E.badatoms(v); ff = E.eqfails(bad)
        print("    APPLIED bad=%s fails=%d SCORE=%d" % (sorted(bad), len(ff), 39033 - len(ff)), flush=True)
        json.dump({"x_%d" % j: int(v[j]) for j in range(E.NV) if v[j] != 0},
                  open('S_reach_%d.json' % (39033 - len(ff)), 'w'))
        print("    WROTE S_reach_%d.json" % (39033 - len(ff)), flush=True)
