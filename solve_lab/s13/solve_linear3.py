#!/usr/bin/env python3
"""
Third pass: put the HANDLES into the unknown set.

solve_linear2.py found the 2-variable system inconsistent at a35760.  That was a
setup artefact, not a fact about the instance:

    a35760 = -x10903*x28961 + x31864 ,  x28961 = p
    a35759 =  5113045*x7075*x9118 - x29854 ,  x29854 = x1329*p

x1329 and x10903 are FREE INPUTS occurring only here -- they are the quotient
handles whose whole job is to absorb (value)/p.  Holding them fixed while asking
whether the check can vanish guarantees a contradiction.

Correct formulation: the unknowns are every FREE INPUT reachable in the cone,
including the handles.  Because handles enter with coefficient p, the honest
condition splits:
    * mod p, the handle terms vanish -> the checks impose congruences on the
      NON-handle unknowns;
    * over Z, each handle then absorbs its check's value / p exactly.
So we solve mod p in the non-handle unknowns, then lift.

Usage: python3 solve_linear3.py
"""
import os, sys, time, json, random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
from reduce_tw import downstream
import fwd_frame as F

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
D2 = [8731, 9118]


def main():
    t0 = time.time()
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    unknown, touched = downstream(D2)
    checks = sorted(a for a in touched if L.atom_out.get(a) is None)
    CLEAN = '--clean' in sys.argv
    if CLEAN:
        frozen = {}          # let every gate hold: the ON-MANIFOLD frame
        bad = []
        print("FRAME: CLEAN (no frozen gates -- all gate atoms must hold)")
    else:
        frozen, bad = F.self_test(v, unknown, touched)
        assert not bad
        print("FRAME: WITNESS (2 gates deliberately broken and frozen)")
    free_all = {t for t in range(L.NVARS) if t not in L.definer}

    # every free input appearing anywhere in the cone's checks / atoms
    occ = defaultdict(int)
    for a in range(len(L.avars)):
        for x in L.avars[a]:
            occ[x] += 1
    cone_free = sorted({x for a in touched for x in L.avars[a] if x in free_all})
    handles = [x for x in cone_free if occ[x] == 1]
    knobs = [x for x in cone_free if occ[x] > 1]
    print(f"cone checks {len(checks)};  free inputs in cone {len(cone_free)}")
    print(f"  HANDLES (occur in exactly one atom, absorb value/p) : "
          f"{len(handles)} {['x%d'%x for x in handles]}")
    print(f"  KNOBS   (genuine unknowns)                          : "
          f"{len(knobs)} {['x%d'%x for x in knobs]}")

    # ---- which checks own a handle => they are absorbable ------------------
    absorbable, hard = [], []
    for a in checks:
        near = set(L.avars[a])
        for x in list(near):
            d = L.definer.get(x)
            if d is not None:
                near |= L.avars[d]
        (absorbable if (near & set(handles)) else hard).append(a)
    print(f"\n  checks with a private handle (absorbable) : {len(absorbable)}")
    print(f"  HARD checks (must vanish on their own)    : {len(hard)}")

    # ---- fit each check as linear mod p in the KNOBS -----------------------
    print("\n" + "=" * 74)
    print("FIT each check as a linear form mod p in the knobs")
    print("=" * 74)
    base = {x: v[x] for x in knobs}

    def ev(delta):
        ov = {x: base[x] + delta.get(x, 0) for x in knobs}
        val = F.evaluate(v, ov, unknown, touched, frozen)
        return {a: F.check_value(a, v, val) for a in checks}

    f0 = ev({})
    cols = {}
    for x in knobs:
        fx = ev({x: 1})
        cols[x] = {a: (fx[a] - f0[a]) % P for a in checks}

    # verify linearity at random points
    random.seed(4)
    goodcnt = defaultdict(int)
    T = 3
    for _ in range(T):
        d = {x: random.randrange(1, 10**9) for x in knobs}
        act = ev(d)
        for a in checks:
            pred = (f0[a] + sum(cols[x][a] * d[x] for x in knobs)) % P
            if pred == act[a] % P:
                goodcnt[a] += 1
    lin = [a for a in checks if goodcnt[a] == T]
    print(f"  checks linear mod p in all knobs : {len(lin)}/{len(checks)}")

    # ---- solve: for absorbable checks require value == 0 mod p ------------
    print("\n" + "=" * 74)
    print("SOLVE the congruence system over the knobs")
    print("=" * 74)
    # NOTE: absorbable checks are NOT dropped.  A handle absorbs (value)/p only
    # when p | value, so the mod-p congruence is still binding.  Excluding them
    # would delete the very constraints that bind.
    rows = []
    for a in lin:
        rows.append(([cols[x][a] % P for x in knobs], (-f0[a]) % P, a))
    print(f"  rows {len(rows)}  x  cols {len(knobs)}")

    # gaussian elimination mod p
    M = [r[0][:] + [r[1], r[2]] for r in rows]
    n = len(knobs)
    piv_rows, rank = [], 0
    for col in range(n):
        pr = None
        for i in range(rank, len(M)):
            if M[i][col] % P != 0:
                pr = i; break
        if pr is None:
            continue
        M[rank], M[pr] = M[pr], M[rank]
        inv = pow(M[rank][col], -1, P)
        M[rank] = [(x * inv) % P for x in M[rank][:n + 1]] + [M[rank][n + 1]]
        for i in range(len(M)):
            if i != rank and M[i][col] % P != 0:
                f = M[i][col] % P
                M[i] = [(M[i][k] - f * M[rank][k]) % P
                        for k in range(n + 1)] + [M[i][n + 1]]
        piv_rows.append(col); rank += 1
    incons = [r for r in M if all(r[k] % P == 0 for k in range(n))
              and r[n] % P != 0]
    print(f"  rank {rank} of {n} columns;  INCONSISTENT rows : {len(incons)}")
    if incons:
        print(f"    from checks {['a%d'%r[n+1] for r in incons[:10]]}")
        print("\n  => the congruence system has NO solution in the knobs of this")
        print("     cone.  This is a linear-algebra fact over GF(p): it is not")
        print("     something a better decomposition, a lower treewidth, or an")
        print("     annealer can change.  To move it one must LEAVE this frame")
        print("     (change which gates are broken / which branch is taken).")
    else:
        sol = {}
        for i, col in enumerate(piv_rows):
            sol[knobs[col]] = M[i][n] % P
        print(f"  CONSISTENT.  solution (mod p) on {len(sol)} knobs")
        json.dump({f'x_{k}': int(t) for k, t in sol.items()},
                  open(os.path.join(HERE, 'knob_solution.json'), 'w'))
        print(f"  written -> {os.path.join(HERE, 'knob_solution.json')}")

    print(f"\n  {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
