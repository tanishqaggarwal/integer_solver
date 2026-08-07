#!/usr/bin/env python3
"""Hunt a compensator for eq8680.

eq8680 is the single row that pins the knob d28730 to zero and so costs the
deliverable its seventh equation.  Its 18 atoms fall into 8 clean PAIRS that each
share a private variable: adding a whole pair to the support turns that shared
variable into a knob, and that knob moves eq8680's core.  So a pair is a candidate
compensator even though neither of its atoms is one on its own.

For each candidate support we compute minfail EXACTLY (enumerate every sacrificed
set within budget), with three sound prunings:
  row == 0 and base == 0  -> equation can never fail, drop it
  row == 0 and base != 0  -> equation can never be satisfied, charge the budget
  otherwise               -> active, enumerate
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


def minfail_exact(E, base, Mat, budget=6, tlimit=600):
    """Smallest number of failing equations, or None if > budget."""
    act = []
    forced = 0
    for i in range(len(E)):
        if all(x == 0 for x in Mat[i]):
            if base[i] != 0:
                forced += 1
        else:
            act.append(i)
    if forced > budget:
        return None, forced, len(act)
    rem = budget - forced
    t0 = time.time()
    for k in range(0, rem + 1):
        for drop in itertools.combinations(act, k):
            keep = [i for i in act if i not in drop]
            if solve_int([Mat[i] for i in keep], [-base[i] for i in keep]) is not None:
                return forced + k, forced, len(act)
            if time.time() - t0 > tlimit:
                return 'timeout', forced, len(act)
    return None, forced, len(act)


def pairs_in(eq):
    ats = [a for _, a in M.eq_terms[eq]]
    share = collections.defaultdict(list)
    for a in ats:
        for x in M.avars[a]:
            share[x].append(a)
    out = []
    for x, lst in share.items():
        if len(lst) == 2 and len(v2a[x]) == 2 and set(v2a[x]) == set(lst):
            out.append((x, tuple(sorted(lst))))
    return sorted(set(out), key=lambda t: t[1])


def main():
    P = pairs_in(EQ)
    print(f"eq{EQ} has {len(M.eq_terms[EQ])} atoms; private-variable PAIRS inside it:")
    coefs = {a: c for c, a in M.eq_terms[EQ]}
    for x, (a1, a2) in P:
        w = list(wit); w[x] = wit[x] + 1
        d1 = M.atom_val(a1, w) - av[a1]
        d2 = M.atom_val(a2, w) - av[a2]
        net = coefs[a1] * d1 + coefs[a2] * d2
        print(f"   knob X{x}: a{a1}(coef {coefs[a1]:+d}, d={d1}) a{a2}(coef {coefs[a2]:+d}, d={d2})"
              f"   net effect on eq{EQ} core = {net}")
    print(flush=True)
    results = []
    t0 = time.time()
    # singles: base + one pair
    for x, pr in P:
        SUP = BASE + list(pr)
        r = build(SUP)
        if r is None:
            print(f"  pair {pr}: non-linear knob, skipped"); continue
        SUPs, knobs, E, base, Mat = r
        mf, forced, nact = minfail_exact(E, base, Mat)
        results.append({'add': list(pr), 'nE': len(E), 'knobs': len(knobs),
                        'minfail': mf, 'forced': forced, 'active': nact})
        tag = f"minfail={mf}" if mf is not None else "minfail > 6"
        star = '   *** BEATS 39,026' if isinstance(mf, int) and mf < 7 else ''
        print(f"  +pair a{pr[0]},a{pr[1]}: |E|={len(E)} knobs={len(knobs)} "
              f"forced_fail={forced} active={nact}  {tag}{star}", flush=True)
    # pairs of pairs
    print("\ntwo pairs at once:", flush=True)
    for (x1, p1), (x2, p2) in itertools.combinations(P, 2):
        SUP = BASE + list(p1) + list(p2)
        r = build(SUP)
        if r is None:
            continue
        SUPs, knobs, E, base, Mat = r
        mf, forced, nact = minfail_exact(E, base, Mat, tlimit=240)
        results.append({'add': list(p1) + list(p2), 'nE': len(E),
                        'knobs': len(knobs), 'minfail': mf, 'forced': forced,
                        'active': nact})
        tag = f"minfail={mf}" if mf is not None else "minfail > 6"
        star = '   *** BEATS 39,026' if isinstance(mf, int) and mf < 7 else ''
        print(f"  +a{p1[0]},a{p1[1]} +a{p2[0]},a{p2[1]}: |E|={len(E)} knobs={len(knobs)} "
              f"forced={forced} active={nact}  {tag}{star}", flush=True)
    json.dump(results, open(os.path.join(HERE, 'eq8680_result.json'), 'w'), indent=1)
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
