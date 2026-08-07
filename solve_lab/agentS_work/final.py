"""Configuration-dependence check: at the configurations where a7389/a10187/a20212/a28647 all
   reach 0 mod p, re-measure knobs+handles+targets and compute the exact reachable lattice on
   a20215.  This is the one row the BFS never drives to 0."""
import sys, json, glob, collections, pickle, time
sys.path.insert(0, '.')
import common as C, lat2
import harness as H, engine as E, fast, sparse, intsolve
P = C.P

hits = sorted(glob.glob('bfs2_hit_*.json')) + sorted(glob.glob('bfs_hit_*.json'))
print("hit configurations:", len(hits), flush=True)
done = 0
for fn in hits:
    a = {int(k): int(v) for k, v in json.load(open(fn)).items()}
    seed = dict(C.BASE)
    seed.update(a)
    v0 = E.forward(seed)
    ns = {'v': v0, '__builtins__': {}}
    rows = {r: eval(H.acodes[r], ns) % P for r in C.ROWS}
    nz = [r for r in C.ROWS if rows[r]]
    if nz != [20215]:
        continue
    done += 1
    print("\n=== %s : ONLY a20215 nonzero mod p (= %s) ===" % (fn, rows[20215]), flush=True)
    t0 = time.time()
    v0, bad0, aff, atoms, hs = lat2.system(seed)
    nh = [x for x in atoms if x not in hs]
    print("   bad=%s  affine knobs=%d  atoms=%d  no-handle=%s (%.0fs)"
          % (sorted(bad0), len(aff), len(atoms), nh, time.time() - t0), flush=True)
    print("   a20215 handle step: %s" % ('p' if hs.get(20215) and abs(hs[20215][1]) == P
                                         else (hs.get(20215) and hs[20215][1])), flush=True)
    knobs = sorted(aff)
    other = [x for x in atoms if x != 20215]
    A = [[aff[f].get(x, 0) for f in knobs] for x in other]
    b = [-bad0.get(x, 0) for x in other]
    n0, ker = intsolve.solve_int(A, b)
    if n0 is None:
        print("   other rows INFEASIBLE (kernel dim %d)" % len(ker), flush=True)
        continue
    Dt = [aff[f].get(20215, 0) for f in knobs]
    r = bad0.get(20215, 0) + sum(Dt[j] * n0[j] for j in range(len(knobs)))
    gens = sorted({abs(sum(Dt[j] * k[j] for j in range(len(knobs)))) for k in ker} - {0})
    import math
    g = 0
    for x in gens:
        g = math.gcd(g, x)
    print("   other rows FEASIBLE (kernel dim %d); reachable step on a20215 = gcd = %s"
          % (len(ker), 'p' if g == P else ('%d*p' % (g // P) if g and g % P == 0 else g)), flush=True)
    print("   a20215 residual at particular solution mod that step: %s"
          % (r % g if g else r), flush=True)
    print("   -> a20215 CAN be zeroed? %s" % (bool(g) and r % g == 0), flush=True)
    if done >= 4:
        break
print("\nchecked %d configurations" % done, flush=True)
