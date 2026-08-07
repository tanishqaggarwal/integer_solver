#!/usr/bin/env python3
"""Cascade (two-level) knobs.

A variable u whose atoms all lie in the support is a knob (already used).  A
variable u with ONE atom b outside the support is ALSO usable if b can be held at
zero by re-solving it for some other variable Y whose own atoms lie inside
support u {b}.  Then b never becomes nonzero, so it costs NO equations, while u
becomes a new knob.  Closing this rule to a fixed point enlarges the realisable
parameter space without enlarging E(S).

This is exactly the move class that a one-level knob analysis cannot see.
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


def unit_solvable(a, y):
    """Is atom a exactly solvable for variable y over Z (coefficient +-1, y not squared)?"""
    q = M.polys[a]
    c1 = 0
    for m, c in q.items():
        if y in m:
            if len(m) != 1:
                return False
            c1 += c
    return abs(c1) == 1


def close(SUP, maxrounds=12):
    """Return (R, dependent, knobs).  R = atoms held at zero by a dependent var."""
    T = set(SUP)
    R = set()
    dep = {}                       # variable -> atom that keeps it pinned
    for _ in range(maxrounds):
        knobs = {x for x in v2a if x not in dep and set(v2a[x]) <= (T | R)}
        added = False
        for b in range(M.na):
            if b in T or b in R:
                continue
            # candidate dependent variables inside b
            for y in M.avars[b]:
                if y in dep or y in knobs:
                    continue
                if set(v2a[y]) - (T | R | {b}):
                    continue
                if not unit_solvable(b, y):
                    continue
                R.add(b); dep[y] = b; added = True
                break
        if not added:
            break
    knobs = sorted({x for x in v2a if x not in dep and set(v2a[x]) <= (T | R)})
    return R, dep, knobs


def resolve(w, R, dep, rounds=40):
    """Re-solve every repairable atom for its dependent variable, over Z."""
    order = list(dep.items())
    for _ in range(rounds):
        changed = False
        for y, b in order:
            q = M.polys[b]
            c1 = 0; rest = 0
            ok = True
            for m, c in q.items():
                if y in m:
                    if len(m) != 1:
                        ok = False; break
                    c1 += c
                else:
                    t = c
                    for x in m:
                        t *= w[x]
                    rest += t
            if not ok or c1 == 0:
                continue
            if rest % c1:
                return False
            nv = -rest // c1
            if w[y] != nv:
                w[y] = nv; changed = True
        if not changed:
            return True
    return True


def main():
    SUP = sorted(RES)
    R, dep, knobs = close(SUP)
    E = sorted({e for a in SUP for e, _ in M.atom_eqs[a]})
    print(f"support {len(SUP)} atoms, |E|={len(E)}")
    print(f"cascade-repairable atoms held at zero: {len(R)}")
    print(f"dependent variables: {len(dep)}")
    print(f"KNOBS after cascade closure: {len(knobs)}  {knobs}", flush=True)
    # verify each knob acts linearly on the support after re-solving
    D = {}
    bad = []
    for k in knobs:
        w = list(wit); w[k] = wit[k] + 1
        if not resolve(w, R, dep):
            bad.append(k); continue
        # every atom outside the support must still be zero
        d1 = [M.atom_val(a, w) - av[a] for a in SUP]
        leak = [b for b in R if M.atom_val(b, w) != 0]
        if leak:
            bad.append(k); continue
        w2 = list(wit); w2[k] = wit[k] + 2
        resolve(w2, R, dep)
        d2 = [M.atom_val(a, w2) - av[a] for a in SUP]
        if any(d2[i] != 2 * d1[i] for i in range(len(SUP))):
            bad.append(k); continue
        D[k] = d1
    print(f"usable linear knobs: {len(D)}   rejected: {len(bad)}", flush=True)
    if not D:
        return
    ks = sorted(D)
    base = []; Mat = []
    for e in E:
        coef = {a: c for c, a in M.eq_terms[e]}
        base.append(sum(coef.get(a, 0) * av[a] for a in SUP))
        Mat.append([sum(coef.get(SUP[i], 0) * D[k][i] for i in range(len(SUP)))
                    for k in ks])
    n = len(E)
    t0 = time.time()
    hit = None
    for kk in range(0, 7):
        for drop in itertools.combinations(range(n), kk):
            keep = [i for i in range(n) if i not in drop]
            if solve_int([Mat[i] for i in keep], [-base[i] for i in keep]) is not None:
                hit = (kk, [E[i] for i in drop]); break
        if hit:
            break
    if hit:
        print(f"*** minfail = {hit[0]}  sacrificing {hit[1]}  => SCORE {M.ne - hit[0]}")
    else:
        print(f"PROVED minfail >= 7 with {len(D)} cascade knobs "
              f"(all C({n},<=6) subsystems integrally unsolvable)")
    print(f"{time.time()-t0:.0f}s")
    json.dump({'nknobs': len(D), 'knobs': ks, 'E': E,
               'minfail': hit[0] if hit else '>=7'},
              open(os.path.join(HERE, 'cascade_result.json'), 'w'))


if __name__ == '__main__':
    main()
