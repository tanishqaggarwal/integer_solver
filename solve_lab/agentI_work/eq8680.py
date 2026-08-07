#!/usr/bin/env python3
"""Hunt a compensator for eq8680 -- the single row that pins the knob d28730 to
zero and so costs the deliverable its seventh equation.

General knob rule: for a variable x, adding ALL of x's atoms to the support makes
x a knob.  So the candidate compensator GROUPS are exactly the sets v2a[x] for
variables x that occur in eq8680's atoms.  A group is admissible if its atoms'
own equations can still be paid for.

minfail is computed EXACTLY by depth-first branch and bound over "include this
equation / sacrifice it", pruning as soon as the included system becomes
integrally unsolvable or the sacrifice budget is exceeded.
"""
import collections, itertools, json, os, sys, time
from model import Model, load_assign
from intsolve import solve_int

HERE = os.path.dirname(os.path.abspath(__file__))
M = Model()
wit = load_assign(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
av = [M.atom_val(a, wit) for a in range(M.na)]
v2a = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        v2a[x].append(i)
RES = [23432, 23433, 36225, 36226, 36227, 36228, 36229, 37537, 37538]
BASE = RES + [23434]
EQ = 8680


def build(SUP):
    SUP = sorted(set(SUP))
    allv = set()
    for a in SUP:
        allv |= M.avars[a]
    knobs = [x for x in sorted(allv) if all(a in SUP for a in v2a[x])]
    if not knobs:
        return None
    D = {}
    for k in knobs:
        w = list(wit); w[k] = wit[k] + 1
        d1 = [M.atom_val(a, w) - av[a] for a in SUP]
        w[k] = wit[k] + 2
        d2 = [M.atom_val(a, w) - av[a] for a in SUP]
        if any(d2[i] != 2 * d1[i] for i in range(len(SUP))):
            return None
        D[k] = d1
    E = sorted({e for a in SUP for e, _ in M.atom_eqs[a]})
    base = []; Mat = []
    for e in E:
        coef = {a: c for c, a in M.eq_terms[e]}
        base.append(sum(coef.get(a, 0) * av[a] for a in SUP))
        Mat.append([sum(coef.get(SUP[i], 0) * D[k][i] for i in range(len(SUP)))
                    for k in knobs])
    return SUP, knobs, E, base, Mat


def minfail_bnb(E, base, Mat, budget=6, tlimit=900):
    """Exact minimum number of failing equations, capped at `budget`."""
    act = []
    forced = 0
    for i in range(len(E)):
        if all(x == 0 for x in Mat[i]):
            if base[i] != 0:
                forced += 1
        else:
            act.append(i)
    if forced > budget:
        return None, forced, len(act), 0
    rem = budget - forced
    # order: rows already satisfied at d=0 first (cheap to keep), then the rest
    act.sort(key=lambda i: (base[i] != 0))
    n = len(act)
    best = [rem + 1]
    nodes = [0]
    t0 = time.time()

    def rec(i, rows, rhs, exc):
        if time.time() - t0 > tlimit:
            raise TimeoutError
        nodes[0] += 1
        if exc >= best[0]:
            return
        if i == n:
            best[0] = exc
            return
        j = act[i]
        r2 = rows + [Mat[j]]; b2 = rhs + [-base[j]]
        if solve_int(r2, b2) is not None:
            rec(i + 1, r2, b2, exc)
        if exc + 1 < best[0]:
            rec(i + 1, rows, rhs, exc + 1)

    try:
        rec(0, [], [], 0)
    except TimeoutError:
        return 'timeout', forced, n, nodes[0]
    if best[0] > rem:
        return None, forced, n, nodes[0]
    return forced + best[0], forced, n, nodes[0]


def groups_for(eq):
    """Candidate compensator groups: v2a[x] for variables x inside eq's atoms,
    restricted to groups not already fully in the base support."""
    ats = set(a for _, a in M.eq_terms[eq])
    out = {}
    for a in ats:
        for x in M.avars[a]:
            g = tuple(sorted(v2a[x]))
            if set(g) <= set(BASE):
                continue
            out.setdefault(g, set()).add(x)
    return out


def report(tag, SUP, budget=6, tlimit=900):
    r = build(SUP)
    if r is None:
        print(f"  {tag}: non-linear knob -> skipped", flush=True)
        return None
    SUPs, knobs, E, base, Mat = r
    mf, forced, nact, nodes = minfail_bnb(E, base, Mat, budget, tlimit)
    if mf is None:
        txt = "minfail > 6"
    elif mf == 'timeout':
        txt = f"TIMEOUT after {nodes} nodes"
    else:
        txt = f"minfail = {mf}"
    star = '   *** BEATS 39,026' if isinstance(mf, int) and mf < 7 else ''
    print(f"  {tag}: |E|={len(E)} knobs={len(knobs)} forced={forced} "
          f"active={nact} nodes={nodes}  {txt}{star}", flush=True)
    return {'tag': tag, 'SUP': SUPs, 'nE': len(E), 'knobs': len(knobs),
            'minfail': mf, 'forced': forced, 'active': nact}


def main():
    G = groups_for(EQ)
    coefs = {a: c for c, a in M.eq_terms[EQ]}
    print(f"eq{EQ}: {len(M.eq_terms[EQ])} atoms.  Candidate knob GROUPS "
          f"(all atoms of one variable):")
    cand = []
    for g, xs in sorted(G.items(), key=lambda kv: len(kv[0])):
        inside = [a for a in g if a in coefs]
        if not inside:
            continue
        x = min(xs)
        w = list(wit); w[x] = wit[x] + 1
        net = sum(coefs[a] * (M.atom_val(a, w) - av[a]) for a in inside)
        newatoms = [a for a in g if a not in BASE]
        cand.append((len(newatoms), net, g, x))
        print(f"   knob X{x}: group {g} ({len(newatoms)} new atoms) "
              f"net effect on eq{EQ} core = {net}")
    cand.sort()
    print(flush=True)
    res = []
    print("single groups added to the deliverable's support + a23434:", flush=True)
    for nnew, net, g, x in cand:
        if net == 0:
            print(f"  group {g}: net effect 0 on eq{EQ} -- cannot compensate", flush=True)
            continue
        r = report(f"+X{x} group {g}", BASE + list(g))
        if r:
            res.append(r)
        json.dump(res, open(os.path.join(HERE, 'eq8680_result.json'), 'w'), indent=1)
    print("\npairs of groups:", flush=True)
    usable = [(nnew, net, g, x) for nnew, net, g, x in cand if net != 0]
    for (n1, t1, g1, x1), (n2, t2, g2, x2) in itertools.combinations(usable, 2):
        r = report(f"+X{x1},X{x2}", BASE + list(g1) + list(g2), tlimit=400)
        if r:
            res.append(r)
        json.dump(res, open(os.path.join(HERE, 'eq8680_result.json'), 'w'), indent=1)
    print("\ndone", flush=True)


if __name__ == '__main__':
    main()

