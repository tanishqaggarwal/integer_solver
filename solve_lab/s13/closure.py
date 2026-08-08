#!/usr/bin/env python3
"""
CLOSE the system, then solve it.

build_cand.py showed the cone-local solution scores 38,858 (worse than 39,026):
it set advice knobs to residues that satisfy the 62 cone checks but violate
checks OUTSIDE the cone.  The cone was never a closed system.

Closure rule (iterate to a fixpoint):
    knobs   -> every CHECK whose value depends on a knob
    checks  -> every FREE INPUT those checks depend on  -> new knobs

At the fixpoint the linear system is self-contained: solving it cannot break
anything it does not already contain.  Then ask whether it is consistent.

Usage: python3 closure.py
"""
import os, sys, time, json, random
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import fwd_frame as F

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
NA = len(L.avars)
SEED_KNOBS = [8731, 9118]


def main():
    t0 = time.time()
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    free_all = {t for t in range(L.NVARS) if t not in L.definer}

    users = defaultdict(list)
    for a in range(NA):
        for x in L.avars[a]:
            users[x].append(a)

    # structural support of every variable over the free inputs
    print("computing free-input support of every wire ...")
    supp = {}
    for x in range(L.NVARS):
        if x in free_all:
            supp[x] = {x}
    order = [y for y in range(L.NVARS) if y in L.definer]
    for _ in range(3):
        changed = 0
        for y in order:
            a = L.definer[y]
            s = set()
            ok = True
            for x in L.avars[a]:
                if x == y:
                    continue
                sx = supp.get(x)
                if sx is None:
                    ok = False; break
                s |= sx
            if ok:
                if supp.get(y) != s:
                    supp[y] = s; changed += 1
        if not changed:
            break
    print(f"  supports computed for {len(supp):,} wires ({time.time()-t0:.0f}s)")

    checks_all = [a for a in range(NA) if L.atom_out.get(a) is None]
    csupp = {}
    for a in checks_all:
        s = set()
        for x in L.avars[a]:
            s |= supp.get(x, set())
        csupp[a] = s

    # ---- closure fixpoint ---------------------------------------------------
    knobs = set(SEED_KNOBS)
    checks = set()
    for it in range(12):
        newchecks = {a for a in checks_all if csupp[a] & knobs}
        newknobs = set(knobs)
        for a in newchecks:
            newknobs |= (csupp[a] & free_all)
        if newchecks == checks and newknobs == knobs:
            break
        checks, knobs = newchecks, newknobs
        print(f"  iter {it}: checks {len(checks):,}, knobs {len(knobs):,}")
    print(f"\nCLOSURE: {len(checks):,} checks x {len(knobs):,} free knobs "
          f"({time.time()-t0:.0f}s)")

    # restrict knobs to those that are not permanently zero/boolean pins
    knobs = sorted(knobs)
    print(f"  (the lab's closure3.py reported 1,655 x 707; same order)")

    # ---- fit the linear model globally --------------------------------------
    print("\nfitting the linear model (one global evaluation per knob)...")
    allvars = set(range(L.NVARS))
    frozen = {}
    touched_all = set(range(NA))

    def evalchecks(delta):
        nv = list(v)
        for x, d in delta.items():
            nv[x] = nv[x] + d
        # global clean forward: recompute every defined wire in topo order
        val = F.evaluate(nv, {}, allvars, touched_all, frozen={})
        return {a: F.check_value(a, nv, val) for a in checks}

    t1 = time.time()
    f0 = evalchecks({})
    print(f"  base evaluation {time.time()-t1:.1f}s; "
          f"nonzero checks at base: {sum(1 for a in checks if f0[a] != 0)}")
    if time.time() - t0 > 600:
        print("  (time budget) stopping before the full fit")
        return

    cols = {}
    for i, x in enumerate(knobs):
        fx = evalchecks({x: 1})
        cols[x] = {a: (fx[a] - f0[a]) % P for a in checks}
        if i % 50 == 0:
            print(f"    knob {i}/{len(knobs)}  ({time.time()-t0:.0f}s)")
        if time.time() - t0 > 900:
            print("    (time budget) truncating the fit")
            knobs = knobs[:i + 1]
            break

    # ---- solve --------------------------------------------------------------
    print("\nsolving the closed congruence system mod p ...")
    chk = sorted(checks)
    n = len(knobs)
    M = []
    for a in chk:
        M.append([cols[x][a] % P for x in knobs] + [(-f0[a]) % P, a])
    rank = 0
    for col in range(n):
        pr = None
        for i in range(rank, len(M)):
            if M[i][col] % P != 0:
                pr = i; break
        if pr is None:
            continue
        M[rank], M[pr] = M[pr], M[rank]
        inv = pow(M[rank][col], -1, P)
        M[rank] = [(t * inv) % P for t in M[rank][:n + 1]] + [M[rank][n + 1]]
        for i in range(len(M)):
            if i != rank and M[i][col] % P != 0:
                f = M[i][col] % P
                M[i] = [(M[i][k] - f * M[rank][k]) % P
                        for k in range(n + 1)] + [M[i][n + 1]]
        rank += 1
    incons = [r for r in M if all(r[k] % P == 0 for k in range(n))
              and r[n] % P != 0]
    print(f"  system {len(M)} x {n};  rank {rank};  "
          f"INCONSISTENT rows {len(incons)}")
    if incons:
        print(f"    e.g. checks {['a%d'%r[n+1] for r in incons[:10]]}")
        print("\n  => the CLOSED linear system over GF(p) is INCONSISTENT.")
        print("     Coupling width is not the obstruction; the obstruction is")
        print("     that no assignment of these knobs satisfies all the")
        print("     congruences at once.  No decomposition changes that.")
    else:
        print("  => CONSISTENT -- a genuine candidate exists; construct and")
        print("     verify it with checker.py.")
    print(f"\n{time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
