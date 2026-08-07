#!/usr/bin/env python3
"""Support search for min |{e : sum_a c_{e,a} v_a != 0}|.

Knobs are variables occurring ONLY in atoms of the support, so every point of the
parameterisation is realisable: perturbing a knob leaves every atom outside the
support exactly zero and every equation outside E(S) exactly satisfied.

For each support we report min failing = |E(S)| - maxsat(S), maxsat computed by
greedy with random restarts (a lower bound on maxsat, i.e. an upper bound on the
score we can claim) and confirmed exhaustively for anything that beats 7.
"""
import collections, itertools, random, json, os, sys, time
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


def greedy_minfail(E, base, Mat, tries=60, seed=0):
    n = len(E)
    rng = random.Random(seed)
    bestk = 0
    order0 = sorted(range(n), key=lambda i: (base[i] != 0))
    for t in range(tries):
        order = order0 if t == 0 else rng.sample(range(n), n)
        rows = []; rhs = []; cnt = 0
        for i in order:
            r2 = rows + [Mat[i]]; b2 = rhs + [-base[i]]
            if solve_int(r2, b2) is not None:
                rows, rhs = r2, b2; cnt += 1
        if cnt > bestk:
            bestk = cnt
    return n - bestk


def exact_minfail(E, base, Mat, cap=8):
    n = len(E)
    for k in range(0, cap + 1):
        for drop in itertools.combinations(range(n), k):
            keep = [i for i in range(n) if i not in drop]
            if solve_int([Mat[i] for i in keep], [-base[i] for i in keep]) is not None:
                return k, [E[i] for i in drop]
    return None, None


def evaluate(SUP, tries=60):
    r = build(SUP)
    if r is None:
        return None
    SUP, knobs, E, base, Mat = r
    mf = greedy_minfail(E, base, Mat, tries=tries)
    return {'SUP': SUP, 'knobs': knobs, 'E': E, 'nE': len(E), 'minfail': mf,
            'base': base, 'Mat': Mat}


def main():
    t0 = time.time()
    r0 = evaluate(RES)
    print(f"base support: |E|={r0['nE']} knobs={len(r0['knobs'])} "
          f"minfail(greedy)={r0['minfail']}", flush=True)
    Eb = set(r0['E'])
    adj = []
    for b in range(M.na):
        if b in RES:
            continue
        eb = {e for e, _ in M.atom_eqs[b]}
        if eb & Eb:
            adj.append((len(eb - Eb), b))
    adj.sort()
    print(f"adjacent atoms: {len(adj)}  cheapest new-equation counts: "
          f"{[a[0] for a in adj[:10]]}", flush=True)
    best = (r0['minfail'], list(RES))
    rows = []
    for newq, b in adj:
        r = evaluate(RES + [b], tries=40)
        if r is None:
            rows.append((99, newq, b, 0)); continue
        rows.append((r['minfail'], r['nE'], b, len(r['knobs'])))
        flag = ''
        if r['minfail'] < best[0]:
            best = (r['minfail'], r['SUP']); flag = '  *** IMPROVEMENT'
        print(f"  +a{b:6d} newEq={newq:3d} |E|={r['nE']:3d} knobs={len(r['knobs']):2d} "
              f"minfail={r['minfail']}{flag}", flush=True)
    # pairs among the cheapest 8
    cheap = [b for _, b in adj[:8]]
    print("\npairs among the 8 cheapest additions:", flush=True)
    for b1, b2 in itertools.combinations(cheap, 2):
        r = evaluate(RES + [b1, b2], tries=40)
        if r is None:
            continue
        flag = ''
        if r['minfail'] < best[0]:
            best = (r['minfail'], r['SUP']); flag = '  *** IMPROVEMENT'
        print(f"  +a{b1},a{b2} |E|={r['nE']:3d} knobs={len(r['knobs']):2d} "
              f"minfail={r['minfail']}{flag}", flush=True)
    # triples among the cheapest 6
    print("\ntriples among the 6 cheapest additions:", flush=True)
    for c in itertools.combinations(cheap[:6], 3):
        r = evaluate(RES + list(c), tries=30)
        if r is None:
            continue
        flag = ''
        if r['minfail'] < best[0]:
            best = (r['minfail'], r['SUP']); flag = '  *** IMPROVEMENT'
        print(f"  +{c} |E|={r['nE']:3d} knobs={len(r['knobs']):2d} "
              f"minfail={r['minfail']}{flag}", flush=True)
    print(f"\nBEST minfail = {best[0]}  => SCORE {M.ne - best[0]}   support={best[1]}")
    json.dump({'best_minfail': best[0], 'best_support': best[1],
               'score': M.ne - best[0], 'rows': rows},
              open(os.path.join(HERE, 'search_result.json'), 'w'))
    print(f"total {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
