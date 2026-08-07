#!/usr/bin/env python3
"""Search over SUPPORTS for a residual placement with fewer than 7 failing equations.

Parameterisation is realisable BY CONSTRUCTION: a knob is a variable that occurs
only in atoms of the support, so perturbing it leaves every atom outside the
support exactly zero and every equation outside E(S) exactly satisfied.

For a support S:
    knobs K(S)   = vars occurring only in S
    E(S)         = equations touched by S
    core_e(d)    = base_e + sum_k Mat[e][k] d_k          (exact, integer, linear)
    failing(S)   = |E(S)| - maxsat(S)
maxsat is computed by greedy with random restarts (lower bound) and, when it
matters, by exact enumeration of the sacrificed set.
"""
import collections, itertools, random, sys, json, os, time
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


def build(SUP):
    SUP = sorted(set(SUP))
    allvars = set()
    for a in SUP:
        allvars |= M.avars[a]
    knobs = [x for x in sorted(allvars) if all(a in SUP for a in v2a[x])]
    if not knobs:
        return None
    D = {}
    for k in knobs:
        w = list(wit); w[k] = wit[k] + 1
        d1 = [M.atom_val(a, w) - av[a] for a in SUP]
        w[k] = wit[k] + 2
        d2 = [M.atom_val(a, w) - av[a] for a in SUP]
        if any(d2[i] != 2 * d1[i] for i in range(len(SUP))):
            return None                      # non-linear knob, skip support
        D[k] = d1
    E = sorted({e for a in SUP for e, _ in M.atom_eqs[a]})
    base = []; Mat = []
    for e in E:
        coef = {a: c for c, a in M.eq_terms[e]}
        base.append(sum(coef.get(a, 0) * av[a] for a in SUP))
        Mat.append([sum(coef.get(SUP[i], 0) * D[k][i] for i in range(len(SUP)))
                    for k in knobs])
    return SUP, knobs, E, base, Mat


def greedy_maxsat(E, base, Mat, tries=150, seed=0):
    n = len(E)
    rng = random.Random(seed)
    best = []
    order0 = sorted(range(n), key=lambda i: (base[i] != 0))
    for t in range(tries):
        order = order0 if t == 0 else rng.sample(range(n), n)
        cur = []
        rows = []; rhs = []
        for i in order:
            r2 = rows + [Mat[i]]; b2 = rhs + [-base[i]]
            if solve_int(r2, b2) is not None:
                rows, rhs = r2, b2
                cur.append(i)
        if len(cur) > len(best):
            best = cur
    return len(best), sorted(E[i] for i in best), best


def exact_can_reach(E, base, Mat, max_fail):
    """Is there T with |E|-|T| <= max_fail and the system on T solvable?"""
    n = len(E)
    for k in range(0, max_fail + 1):
        for drop in itertools.combinations(range(n), k):
            keep = [i for i in range(n) if i not in drop]
            rows = [Mat[i] for i in keep]; rhs = [-base[i] for i in keep]
            if solve_int(rows, rhs) is not None:
                return True, [E[i] for i in drop]
    return False, None


def evaluate(SUP, tries=150):
    r = build(SUP)
    if r is None:
        return None
    SUP, knobs, E, base, Mat = r
    ms, sat, _ = greedy_maxsat(E, base, Mat, tries=tries)
    return {'SUP': SUP, 'knobs': knobs, 'E': E, 'nE': len(E),
            'maxsat': ms, 'failing': len(E) - ms, 'sat': sat,
            'base': base, 'Mat': Mat}


def main():
    t0 = time.time()
    r0 = evaluate(RES)
    print(f"base support (9 residual atoms): |E|={r0['nE']} maxsat={r0['maxsat']} "
          f"failing={r0['failing']} knobs={r0['knobs']}", flush=True)
    best = r0
    # candidate atoms: every atom appearing in an equation of E(S)
    cand = sorted({a for e in r0['E'] for _, a in M.eq_terms[e]} - set(RES))
    print(f"candidate atoms to add: {len(cand)}", flush=True)
    rows = []
    for k, b in enumerate(cand):
        r = evaluate(RES + [b], tries=60)
        if r is None:
            continue
        rows.append((r['failing'], r['nE'], r['maxsat'], b, len(r['knobs'])))
        if r['failing'] < best['failing']:
            best = r
            print(f"  *** IMPROVEMENT adding a{b}: |E|={r['nE']} maxsat={r['maxsat']} "
                  f"failing={r['failing']}", flush=True)
        if k % 5 == 0:
            print(f"  {k}/{len(cand)} best failing={best['failing']} "
                  f"t={time.time()-t0:.0f}s", flush=True)
    rows.sort()
    print("\ntop single additions (failing, |E|, maxsat, atom, #knobs):")
    for r in rows[:15]:
        print("   ", r)
    json.dump({'base_failing': r0['failing'],
               'best_failing': best['failing'],
               'best_SUP': best['SUP'],
               'top': rows[:40]},
              open(os.path.join(HERE, 'beam_result.json'), 'w'))
    print(f"\nBEST failing = {best['failing']}  => SCORE {M.ne - best['failing']}")


if __name__ == '__main__':
    main()
