"""Apply the exact-lattice + sacrifice machinery in the DELIVERABLE's basin.

The deliverable's residual family (a23616/a23617/a36659..a36664) has much smaller equation
footprints than the cfg0 cluster, so it is the only basin where beating 7 is arithmetically
possible.  Everything is re-measured in this basin: knobs, handles AND targets.
"""
import sys, json, collections, pickle, time, itertools
sys.path.insert(0, '.')
import common as C, lattice as L
import harness as H, engine as E, fast, sparse
P = C.P

FOOT = collections.defaultdict(set)
for e, (issq, outer, terms) in enumerate(H.eqt):
    for c, a in terms:
        if a >= 0:
            FOOT[a].add(e)
NF = {a: len(s) for a, s in FOOT.items()}

d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vD = [0] * E.NV
for k, val in d.items():
    vD[int(k.split('_')[1])] = int(val)
seed = {f: vD[f] for f in sorted(E.FREE)}
v0 = E.forward(seed)
bad0 = E.badatoms(v0)
ff0 = E.eqfails(bad0)
print("basin state: bad=%s (nf %s) fails=%d SCORE=%d"
      % (sorted(bad0), [NF[a] for a in sorted(bad0)], len(ff0), 39033 - len(ff0)), flush=True)

cone = set()
for a in bad0:
    cone |= set(E.cone(a)[1])
cone = sorted(cone)
print("union cone free vars:", len(cone), flush=True)

T = L.delta_table(v0, bad0, cone)
aff = {f: dd for f, (dd, afq) in T.items() if afq and dd and not C.isbool(f)}
atoms = set(bad0)
for dd in aff.values():
    atoms |= set(dd)
cache = {}
hs = L.handles_for(v0, bad0, sorted(atoms), cache)
for a, (f, s) in hs.items():
    if f not in aff:
        aff[f] = {a: s}
atoms = set(bad0)
for dd in aff.values():
    atoms |= set(dd)
atoms = sorted(atoms)
nh = [a for a in atoms if a not in hs]
print("affine knobs=%d atoms=%d no-handle=%s (nf %s)"
      % (len(aff), len(atoms), nh, [NF.get(a) for a in nh]), flush=True)

knobs = sorted(aff)
rowd = {a: {f: aff[f][a] for f in knobs if a in aff[f]} for a in atoms}
rhs = {a: -bad0.get(a, 0) for a in atoms}


def solve(keep):
    return sparse.solve_sparse([rowd[a] for a in keep], [rhs[a] for a in keep],
                               verbose=False, maxcore=800, maxcorebits=400_000)[0]


def apply(sol, tag):
    ns = dict(seed)
    for f, dv in sol.items():
        if dv:
            ns[f] = v0[f] + dv
    v = E.forward(ns)
    bad = E.badatoms(v)
    ff = E.eqfails(bad)
    print("  %s -> bad=%s (nf %s) fails=%d SCORE=%d"
          % (tag, sorted(bad), [NF.get(a) for a in sorted(bad)], len(ff), 39033 - len(ff)), flush=True)
    if len(ff) < 7:
        json.dump({"x_%d" % j: int(v[j]) for j in range(E.NV) if v[j] != 0},
                  open('S_basin5_%d.json' % (39033 - len(ff)), 'w'))
        json.dump({str(x): str(int(y)) for x, y in ns.items()}, open('S_basin5_seed.json', 'w'))
        print("  *** WROTE S_basin5_%d.json ***" % (39033 - len(ff)), flush=True)
    return len(ff)


sol = solve(atoms)
print("FULL SYSTEM:", "FEASIBLE!!!" if sol is not None else "infeasible", flush=True)
if sol is not None:
    apply(sol, 'full')
else:
    cand = sorted(atoms, key=lambda a: (NF.get(a, 99), a))
    print("lowest-footprint atoms in play:", [(NF.get(a), a) for a in cand[:16]], flush=True)
    best = 10 ** 9
    for k in (4, 5):
        for S in itertools.combinations(cand[:13], k):
            s2 = solve([a for a in atoms if a not in S])
            if s2 is None:
                continue
            n = apply(s2, 'drop%s nf=%s' % (list(S), [NF.get(a) for a in S]))
            best = min(best, n)
        if best <= 6:
            break
    print("BEST", 39033 - best, flush=True)
